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


def load_chemicals_from_parquet(parquet_path)
    # ~ parquet_path: str,
# ~ ) -> dict[str, tuple[str | None, str | None]]:
    """Read the Parquet file and build {name -> (chebi_iri, pubchem_iri)}.

    The file produced by ChemProp (or similar) is expected to contain at least
    the columns 'name' and 'cid'.  If CHEBI identifiers are unavailable, we set
    them to ``None`` and rely on PubChem mapping only.
    """
    df = pd.read_parquet(parquet_path)

    # ~ if "name" not in df.columns or "cid" not in df.columns:
        # ~ raise ValueError("Parquet file must contain 'name' and 'cid' columns.")

    # ~ # Drop duplicates to avoid redundant SPARQL/HDT look‑ups
    # ~ df_unique = df.drop_duplicates(subset=["name", "cid"])

    # ~ mapping: dict[str, tuple[str | None, str | None]] = {}
    # ~ for _, row in df_unique.iterrows():
        # ~ name: str = str(row["name"]).strip()
        # ~ cid = row["cid"]
        # ~ pubchem_iri = (
            # ~ f"https://identifiers.org/pubchem.compound/{int(cid)}" if pd.notna(cid) else None
        # ~ )
        # ~ mapping[name] = (None, pubchem_iri)  # No CHEBI IRI available
    # ~ return mapping

    chemicals = df.inchi.unique()


def map_chemicals_to_pathways(
    chemicals: dict[str, tuple[str | None, str | None]],
    hdt,
) -> dict[str, set[str]]:
    """Return a mapping {chemical_name -> set(pathway_iri)}.

    Both CHEBI and PubChem identifiers are tried when available.
    """
    chem_to_path: dict[str, set[str]] = collections.defaultdict(set)

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
    df = get_gene_products(pathway, store)
    store.close()
    print(df)

    # ~ # Identify the relevant pathways by chemical
    # ~ chem_to_path = map_chemicals_to_pathways(chemicals, store)
    # ~ unique_pathways = {pwy for pwys in chem_to_path.values() for pwy in pwys}

    # ~ print(unique_pathways)

    # ~ # Get gene products from these pathways
    # ~ path_to_genes = map_pathways_to_genes(unique_pathways, hdt)

    # ~ # Load the property file
    # ~ fname = "predicted_property_names.txt"
    # ~ with open(fname, 'r') as file:
        # ~ properties = file.readlines()
    
    # ~ # Use LLM to predict the most relevant properties from the list
    # ~ response = client.models.generate_content(
        # ~ model = "gemini-2.0-flash",
        # ~ contents = [
            # ~ "Analyze the relevant properties of these pathways:",
            # ~ unique_pathways,
            # ~ "given this list of possible properties:",
            # ~ properties,
        # ~ ],
    # ~ )


if __name__ == "__main__":
    main()
