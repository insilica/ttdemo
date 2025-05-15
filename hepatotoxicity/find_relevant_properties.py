# TODO: use cache
# TODO: split into different files
# TODO: pathway visualization
from tqdm import tqdm
import pandas as pd
import biobricks as bb
from rdflib import Graph, Namespace
from rdflib.plugins.stores import sparqlstore
from rdflib_hdt import HDTStore
import google
from google import genai
import time

import os
from dotenv import load_dotenv
import json, re

from pydantic import BaseModel

class GenePropSchema(BaseModel):
    property_names : list[str]

# import pathlib
# import toxindex.utils.chemprop as chemprop
# import toxindex.utils.simplecache as simplecache

# # -- RDF PREDICATES -----------------------------------------
# CHEBI_PRED   = "http://vocabularies.wikipathways.org/wp#bdbChEBI"
# PUBCHEM_PRED = "http://vocabularies.wikipathways.org/wp#bdbPubChem"
# ISPART_PRED  = "http://purl.org/dc/terms/isPartOf"

# ENTREZ_PRED  = "http://vocabularies.wikipathways.org/wp#bdbEntrezGene"
# UNIPROT_PRED = "http://vocabularies.wikipathways.org/wp#bdbUniprot"
# LABEL_PRED   = "http://www.w3.org/2000/01/rdf-schema#label"


# def map_chemicals_to_pathways(
#     chemicals: dict[str, tuple[str | None, str | None]],
#     hdt,
# ) -> dict[str, set[str]]:
#     """Return a mapping {chemical_name -> set(pathway_iri)}.

#     Both CHEBI and PubChem identifiers are tried when available.
#     """
#     chem_to_path: dict[str, set[str]] = collections.defaultdict(set)

#     for name, (chebi_iri, pubchem_iri) in chemicals.items():
#         id_pairs = (
#             (CHEBI_PRED, chebi_iri),
#             (PUBCHEM_PRED, pubchem_iri),
#         )
#         for pred, iri in id_pairs:
#             if not iri:
#                 continue  # Identifier missing → skip
#             for node, _, _ in hdt.search_triples("", pred, iri):
#                 for _, _, pwy in hdt.search_triples(node, ISPART_PRED, "")[0]:
#                     chem_to_path[name].add(pwy)
#     return chem_to_path


# def map_pathways_to_genes(
#     pathways: set[str],
#     hdt,
# ) -> dict[str, list[str]]:
#     """Return {pathway_iri -> sorted list of human gene labels}."""
#     path_to_genes: dict[str, list[str]] = {}

#     for pwy in pathways:
#         genes: set[str] = set()
#         for gene_pred in (ENTREZ_PRED, UNIPROT_PRED):
#             for node, _, _ in hdt.search_triples("", gene_pred, ""):
#                 if hdt.triple_exist(node, ISPART_PRED, pwy):
#                     for _, _, label in hdt.search_triples(node, LABEL_PRED, "")[0]:
#                         genes.add(label)
#         path_to_genes[pwy] = sorted(genes)
#     return path_to_genes


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
    df = pd.DataFrame(results, columns=["geneProduct", "geneLabel"])
    return df


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


def main():
    # Load environment variables from .env file in the current directory
    load_dotenv()
    client = genai.Client(api_key = os.environ["GEMINI_API_KEY"])

    # # Process the parquet of hepatotoxic chemicals
    # chemicals = load_chemicals_from_parquet("chemprop_predictions.parquet")
    # print(chemicals)

    # # Predict properties and process the predictions
    # cachedir = pathlib.Path('cache') / 'function_cache' / 'chemprop_predictions_cache'
    # cachedir.mkdir(parents=True, exist_ok=True)
    # pred_func = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)
    # for chemical in chemicals:
    #     prediction = pred_func(chemical)

    # Get all predictions for all chemicals included in parquet
    predictions = pd.read_parquet("chemprop_predictions.parquet")
    # # check if all property_title entries are strings
    # if not predictions['property_title'].apply(lambda x: isinstance(x, str)).all():
    #     raise ValueError("Not all property_title entries are strings.")

    # Load wikipathways biobrick
    wikipathways = bb.assets('wikipathways')
    hdt_path = wikipathways.wikipathways_hdt
    store = HDTStore(hdt_path)    

    # Query a pathway to get the gene products
    pathway = "WP3657"  # Chosen pathway
    pathway_name = " Hematopoietic stem cell gene regulation by GABP alpha/beta complex"
    df = get_gene_products(pathway, store)
    store.close()

    # Load the property file
    fname = "predicted_property_names.txt"
    with open(fname, 'r') as file:
        properties = file.readlines()

    # Remove newline characters and strip whitespace
    properties = [prop.strip() for prop in properties]
    rows = []
    for gene in tqdm(df['geneLabel'], desc = "Predicting properties for genes", unit = "gene"):
        # Use LLM to predict the most relevant properties from the list
        while True:
            try:
                response = client.models.generate_content(
                    model = "gemini-2.0-flash",
                    contents = [
                        # f"Gene: {gene}",
                        # f"Pathway: {pathway_name}",
                        # (
                        #     "From the list below, return a **JSON object** with a single key "
                        #     '"property_names". Its value must be a JSON array containing ONLY '
                        #     "those property strings that are relevant. Use the property names "
                        #     "exactly as given—verbatim—and do not add keys or explanations."
                        # ),
                        # "List of possible properties:",
                        # properties,
                        f"select properties from the list below that are strongly associated with the gene product {gene} in the biological pathway {pathway}\n\n{properties}\n\nOutput one unique property per line and nothing else; copy the property names exactly."
                    ],
                    config = {
                        # "response_mime_type": "application/json",
                        # "response_schema": GenePropSchema,
                        "temperature" : 0.0,  # deterministic output
                    }
                )
                break  # exit the loop if successful
            except google.genai.errors.ServerError as e:
                print(f"Server error: {e}. Retrying...")
                time.sleep(10)  # wait before retrying
        
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
    gene_property_predictions.to_parquet("gene_property_predictions.parquet")

    # select the chemicals with the relevant properties
    rows = []
    # for gene_pathway, property_list in tqdm(gene_property_predictions.items(), desc = "Selecting chemicals with relevant properties"):
    for gpp in tqdm(gene_property_predictions.iterrows(), desc = "Selecting chemicals with relevant properties"):
        # gene, pathway = gene_pathway
        chemicals_values = []
        for pred in predictions.iterrows():
            if pred['property_title'] in gpp['property_list']:
                chemicals_values.append((pred['inchi'], pred['value']))
        # add the gene and its chemicals to the dataframe
        rows.append({
            'gene': gpp['gene'],
            'pathway': gpp['pathway'],
            'chemicals_values': chemicals_values
        })
    gene_property_chemicals = pd.DataFrame(rows)
    # save the dataframe to a parquet file
    gene_property_chemicals.to_parquet("gene_property_chemicals.parquet")

    # group by genes and sum over values
    merged_df = df.groupby("gene", as_index=False).agg({
        "chemicals_values": lambda lsts: sum(lsts, [])  # flatten list of lists
    })

    def compute_sum_and_mean(chem_vals):
        values = [val for _, val in chem_vals]
        return pd.Series({
            "sum_value": sum(values),
            "mean_value_over_chemicals": sum(values) / len(values) if values else float('nan')
        })

    summary_df = merged_df.copy()
    summary_df[["sum_value", "mean_value_over_chemicals"]] = merged_df["chemicals_values"].apply(compute_sum_and_mean)

    # save the summary dataframe to a parquet file
    summary_df.to_parquet("gene_property_chemicals_summary.parquet")
    # print the summary dataframe
    print(summary_df)


if __name__ == "__main__":
    main()
