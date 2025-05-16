from tqdm import tqdm
import pandas as pd
import biobricks as bb
from rdflib import Graph
from rdflib.plugins.stores import sparqlstore
from rdflib_hdt import HDTStore
import google
from google import genai
from tenacity import retry, wait_fixed, retry_if_exception_type
import time
import argparse
import numpy as np

import os
from dotenv import load_dotenv
import json, re

from pydantic import BaseModel

import pathlib

class GenePropSchema(BaseModel):
    property_names : list[str]

def get_gene_products(
    pathway: str,
    store: HDTStore,
):
    g = Graph(store = store)
    query = f"""
    PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>
    PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX dcterms: <http://purl.org/dc/terms/>

    SELECT DISTINCT ?geneProduct ?geneLabel
    WHERE {{
    # --- Pathway node ----------------------------------------------------------
    ?pathway a wp:Pathway ;
            dcterms:identifier "{pathway}" .

    # --- Gene products that belong to that pathway -----------------------------
    ?geneProduct a wp:GeneProduct ;
                dcterms:isPartOf ?pathway ;
                rdfs:label      ?geneLabel .
    }}
    ORDER BY ?geneLabel
    """
    results = g.query(query)
    gene_products = pd.DataFrame(results, columns=["geneProduct", "geneLabel"])
    return gene_products


def get_response_json(response):
    """Alter the response to be a valid JSON object."""
    if isinstance(response, str):
        raw = response
    else:
        raw = response.text

    try:
        data = json.loads(raw)  # the normal path
    except json.JSONDecodeError:
        # quick fixer for the common case of a missing closing brace/bracket
        # 1) keep only the JSON‑looking part
        snippet = re.search(r'\{.*', raw, re.S).group(0)

        # 2) add missing quotes around keys
        if snippet.count('"') % 2 == 1:
            snippet += '"'

        # 3) add a missing closing brace/bracket if obvious
        if snippet.count('[') > snippet.count(']'):
            snippet += ']'
        if snippet.count('{') > snippet.count('}'):
            snippet += '}'

        # 4) try again
        data = json.loads(snippet)

    return data


@ retry(
    retry = retry_if_exception_type(google.genai.errors.ServerError),
    wait = wait_fixed(10),
    reraise = True
)
def get_property_response(client, gene, pathway, properties):
    """Retrying Gemini API call on server error."""
    response = client.models.generate_content(
        model = "gemini-2.0-flash",
        contents = [
            f"select properties from the list below that are strongly associated with the gene product {gene} in the biological pathway {pathway}\n\n{properties}\n\nOutput one unique property per line and nothing else; copy the property names exactly."
        ],
        config = {
            "temperature": 0.0,
        }
    )
    return response


def predict_properties_for_genes(outdir, client, pathway, properties, gene_products, use_cache = True):
    """Use the Gemini LLM to predict properties for genes in a pathway.

    The properties are selected from a list of possible properties.
    """
    # check if the results are cached
    fname = outdir / f"{pathway}.gene_property_predictions.parquet"
    if os.path.exists(fname) and use_cache:
        return pd.read_parquet(fname)

    rows = []
    for gene in tqdm(gene_products['geneLabel'], desc = "Predicting properties for genes", unit = "gene"):
        # Use LLM to predict the most relevant properties from the list
        try:
            response = get_property_response(client, gene, pathway, properties)
        except Exception as e:
            print(f"Failed to get response for gene {gene}: {e}")
            continue  # skip this gene
        
        # convert the output to a list of strings
        property_list = response.text.splitlines()
        
        # remove all of the properties not in the list
        # TODO: implement a fuzzy comparison to be less restrictive (e.g., rapidfuzz)
        n_removed = 0
        for prop in property_list:
            if prop not in properties:
                property_list.remove(prop)
                n_removed += 1
        # print(f"Removed {n_removed} properties not in the list. {len(property_list)} remaining properties.")
        # add the gene and its properties to the dictionary
        rows.append({
            'gene': gene,
            'pathway': pathway,
            'property_list': property_list
        })
        # [(gene, pathway)] = property_list

    # save the responses to a parquet file
    gene_property_predictions = pd.DataFrame(rows)
    gene_property_predictions.to_parquet(fname)

    return gene_property_predictions


def gene_chemical_strengths(outdir, pathway, gene_property_predictions, predictions, use_cache = True):
    """Select the relevant chemicals from the predictions.

    The chemicals are selected based on the properties predicted for the genes.
    """
    fname = outdir / f"{pathway}.gene_property_chemicals.parquet"
    if os.path.exists(fname) and use_cache:
        gene_property_chemicals = pd.read_parquet(fname)
    else:
        rows = []
        # for gene_pathway, property_list in tqdm(gene_property_predictions.items(), desc = "Selecting chemicals with relevant properties"):
        for _, gpp in tqdm(gene_property_predictions.iterrows(), desc = "Selecting chemicals with relevant properties", total = len(gene_property_predictions)):
            # gene, pathway = gene_pathway
            chemicals_values = []
            property_list = gpp['property_list'].tolist()
            for _, pred in predictions.iterrows():
                if pred['property_title'] in property_list:
                    chemicals_values.append((pred['inchi'], pred['value']))

            # add the gene and its chemicals to the dataframe
            rows.append({
                'gene': gpp['gene'],
                'pathway': gpp['pathway'],
                'chemicals_values': chemicals_values
            })

        # Serialize list of tuples to JSON
        gene_property_chemicals = pd.DataFrame(rows)
        gene_property_chemicals['chemicals_values'] = gene_property_chemicals['chemicals_values'].apply(json.dumps)
        # save the dataframe to a parquet file
        gene_property_chemicals.to_parquet(fname)

    # return to list of tuples
    gene_property_chemicals['chemicals_values'] = gene_property_chemicals['chemicals_values'].apply(json.loads)
    return gene_property_chemicals


def main(cachedir, projectdir, outdir, pathway, use_cache_predictions = True, use_cache_chemicals = True):
    # Load environment variables from .env file in the current directory
    load_dotenv()
    client = genai.Client(api_key = os.environ["GEMINI_API_KEY"])

    # Get all predictions for all chemicals included in parquet
    predictions = pd.read_parquet(projectdir / "predict_chemicals" / "chemprop_predictions.parquet")

    # Load wikipathways biobrick
    wikipathways = bb.assets('wikipathways')
    hdt_path = wikipathways.wikipathways_hdt
    store = HDTStore(hdt_path)    

    # Query a pathway to get the gene products
    gene_products = get_gene_products(pathway, store)
    store.close()

    # Load the property file
    fname = cachedir / "predicted_property_names.txt"
    with open(fname, 'r') as file:
        properties = file.readlines()

    # Remove newline characters and strip whitespace
    properties = [prop.strip() for prop in properties]

    # Get predictions for the gene products
    gene_property_predictions = predict_properties_for_genes(
        outdir, client, pathway, properties, gene_products, use_cache_predictions                       
    )
    
    # select the chemicals with the relevant properties
    gene_property_chemicals = gene_chemical_strengths(
        outdir, pathway, gene_property_predictions, predictions, use_cache_chemicals
    )

    # group by genes and sum over values
    merged_gene_products = gene_property_chemicals.groupby("gene", as_index=False).agg({
        "chemicals_values": lambda lsts: sum(lsts, [])  # flatten list of lists
    })

    def compute_stats(chem_vals):
        values = [val for _, val in chem_vals]
        return pd.Series({
            "sum_value": sum(values),
            "mean_value": sum(values) / len(values) if values else float('nan'),
            "2-norm": np.linalg.norm(values),
            "3-norm": np.linalg.norm(values, ord = 3),
            "4-norm": np.linalg.norm(values, ord = 4),
            "max_value": max(values) if values else float('nan'),
        })

    summary_gene_products = merged_gene_products.copy()
    summary_gene_products[["sum_value", "mean_value", "2-norm", "3-norm", "4-norm", "max_value"]] = merged_gene_products["chemicals_values"].apply(compute_stats)
    summary_gene_products.drop(columns = ["chemicals_values"], inplace = True)
    print(summary_gene_products)
    # save the summary dataframe to a parquet file
    summary_gene_products.to_parquet(outdir / f"{pathway}.gene_property_chemicals_summary.parquet")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find relevant properties for genes in a pathway.")
    parser.add_argument(
        "--project",
        type = str,
        default = "hepatotoxic",
        help = "Project name",
    )
    parser.add_argument(
        "--pathway",
        type = str,
        required = True,
        help = "Pathway ID to query",
    )
    parser.add_argument(
        "--use_cache_predictions",
        type = bool,
        default = True,
        help = "Use cached predictions for the genes in the pathway",
    )
    parser.add_argument(
        "--use_cache_chemicals",
        type = bool,
        default = True,
        help = "Use cached predictions for the genes in the pathway",
    )
    args = parser.parse_args()

    cachedir = pathlib.Path("cache")
    projectdir = cachedir / "projects" / args.project
    outdir = projectdir / "gene_property_predictions"
    outdir.mkdir(parents = True, exist_ok = True)
    print(f"cachedir = {cachedir}, projectdir = {projectdir}, outdir = {outdir}")

    main(cachedir, projectdir, outdir, args.pathway, args.use_cache_predictions, args.use_cache_chemicals)
