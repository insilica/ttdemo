import collections
import pandas as pd
import biobricks as bb
from rdflib import Graph, Namespace
from rdflib.plugins.stores import sparqlstore
from rdflib_hdt import HDTStore
from google import genai

import os
from dotenv import load_dotenv

import pathlib
import toxindex.utils.chemprop as chemprop
import toxindex.utils.simplecache as simplecache

# -- RDF PREDICATES -----------------------------------------
CHEBI_PRED   = "http://vocabularies.wikipathways.org/wp#bdbChEBI"
PUBCHEM_PRED = "http://vocabularies.wikipathways.org/wp#bdbPubChem"
ISPART_PRED  = "http://purl.org/dc/terms/isPartOf"

ENTREZ_PRED  = "http://vocabularies.wikipathways.org/wp#bdbEntrezGene"
UNIPROT_PRED = "http://vocabularies.wikipathways.org/wp#bdbUniprot"
LABEL_PRED   = "http://www.w3.org/2000/01/rdf-schema#label"


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

    for name, (chebi_iri, pubchem_iri) in chemicals.items():
        id_pairs = (
            (CHEBI_PRED, chebi_iri),
            (PUBCHEM_PRED, pubchem_iri),
        )
        for pred, iri in id_pairs:
            if not iri:
                continue  # Identifier missing → skip
            for node, _, _ in hdt.search_triples("", pred, iri):
                for _, _, pwy in hdt.search_triples(node, ISPART_PRED, "")[0]:
                    chem_to_path[name].add(pwy)
    return chem_to_path


def map_pathways_to_genes(
    pathways: set[str],
    hdt,
) -> dict[str, list[str]]:
    """Return {pathway_iri -> sorted list of human gene labels}."""
    path_to_genes: dict[str, list[str]] = {}

    for pwy in pathways:
        genes: set[str] = set()
        for gene_pred in (ENTREZ_PRED, UNIPROT_PRED):
            for node, _, _ in hdt.search_triples("", gene_pred, ""):
                if hdt.triple_exist(node, ISPART_PRED, pwy):
                    for _, _, label in hdt.search_triples(node, LABEL_PRED, "")[0]:
                        genes.add(label)
        path_to_genes[pwy] = sorted(genes)
    return path_to_genes


def generate_llm_input(inchi : str, pred_func):
    


def main():
    # Load environment variables from .env file in the current directory
    load_dotenv()
    client = genai.Client(api_key = os.environ["GEMINI_API_KEY"])

    # Load wikipathways biobrick
    wikipathways = bb.assets('wikipathways')
    hdt_path = wikipathways.wikipathways_hdt
    store = HDTStore(hdt_path)

    # Process the parquet of hepatotoxic chemicals
    chemicals = load_chemicals_from_parquet("chemprop_predictions.parquet")
    print(chemicals)

    # Predict properties and process the predictions
    cachedir = pathlib.Path('cache') / 'function_cache' / 'chemprop_predictions_cache'
    cachedir.mkdir(parents=True, exist_ok=True)
    pred_func = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)
    for chemical in chemicals:
        # make the prediction
        prediction = pred_func(chemical)

        # get categories and strengths
        categories = []
        strengths = []
        for i in range(len(prediction)):
            prediction_categories = prediction[i]['property']['categories']
            for j in range(len(prediction_categories)):
                categories.append(prediction_categories[j]['category'])
                strengths.append(prediction_categories[j]['strength'])        
        
        # get statistics
        categories = pd.Series(categories)
        strengths = pd.Series(strengths)
        numerical_predictions = pd.DataFrame({'category': categories, 'strength': strengths})

        cat_group = numerical_predictions.groupby('category')['strength']
        means = cat_group.mean()
        stds  = cat_group.std()
        stats = pd.DataFrame({'mean': means, 'std': stds})
        stats = stats.round(2)  # high precision not needed for LLM

        # call the LLM to (make a prediction)
        llm_input = stats.reset_index().to_dict(orient='records')
        

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

import collections
import pandas as pd
import biobricks as bb
from rdflib import Graph, Namespace
from rdflib.plugins.stores import sparqlstore
from rdflib_hdt import HDTStore
from google import genai

import os
from dotenv import load_dotenv

import pathlib
import toxindex.utils.chemprop as chemprop
import toxindex.utils.simplecache as simplecache

# -- RDF PREDICATES -----------------------------------------
CHEBI_PRED   = "http://vocabularies.wikipathways.org/wp#bdbChEBI"
PUBCHEM_PRED = "http://vocabularies.wikipathways.org/wp#bdbPubChem"
ISPART_PRED  = "http://purl.org/dc/terms/isPartOf"

ENTREZ_PRED  = "http://vocabularies.wikipathways.org/wp#bdbEntrezGene"
UNIPROT_PRED = "http://vocabularies.wikipathways.org/wp#bdbUniprot"
LABEL_PRED   = "http://www.w3.org/2000/01/rdf-schema#label"


def load_chemicals_from_parquet(parquet_path):
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

    for name, (chebi_iri, pubchem_iri) in chemicals.items():
        id_pairs = (
            (CHEBI_PRED, chebi_iri),
            (PUBCHEM_PRED, pubchem_iri),
        )
        for pred, iri in id_pairs:
            if not iri:
                continue  # Identifier missing → skip
            for node, _, _ in hdt.search_triples("", pred, iri):
                for _, _, pwy in hdt.search_triples(node, ISPART_PRED, "")[0]:
                    chem_to_path[name].add(pwy)
    return chem_to_path


def map_pathways_to_genes(
    pathways: set[str],
    hdt,
) -> dict[str, list[str]]:
    """Return {pathway_iri -> sorted list of human gene labels}."""
    path_to_genes: dict[str, list[str]] = {}

    for pwy in pathways:
        genes: set[str] = set()
        for gene_pred in (ENTREZ_PRED, UNIPROT_PRED):
            for node, _, _ in hdt.search_triples("", gene_pred, ""):
                if hdt.triple_exist(node, ISPART_PRED, pwy):
                    for _, _, label in hdt.search_triples(node, LABEL_PRED, "")[0]:
                        genes.add(label)
        path_to_genes[pwy] = sorted(genes)
    return path_to_genes


def generate_llm_input(inchi : str, pred_func):
    # make the prediction
    prediction = pred_func(chemical)

    # get categories and strengths
    categories = []
    strengths = []
    for i in range(len(prediction)):
        prediction_categories = prediction[i]['property']['categories']
        for j in range(len(prediction_categories)):
            categories.append(prediction_categories[j]['category'])
            strengths.append(prediction_categories[j]['strength'])        
    
    # get statistics
    categories = pd.Series(categories)
    strengths = pd.Series(strengths)
    numerical_predictions = pd.DataFrame({'category': categories, 'strength': strengths})

    cat_group = numerical_predictions.groupby('category')['strength']
    means = cat_group.mean()
    stds  = cat_group.std()
    stats = pd.DataFrame({'mean': means, 'std': stds})
    stats = stats.round(2)  # high precision not needed for LLM

    # call the LLM to (make a prediction)
    llm_input = stats.reset_index().to_dict(orient='records')

    return llm_input


def main():
    # Load environment variables from .env file in the current directory
    load_dotenv()
    client = genai.Client(api_key = os.environ["GEMINI_API_KEY"])

    # Load wikipathways biobrick
    wikipathways = bb.assets('wikipathways')
    hdt_path = wikipathways.wikipathways_hdt
    store = HDTStore(hdt_path)

    # Process the parquet of hepatotoxic chemicals
    chemicals = load_chemicals_from_parquet("chemprop_predictions.parquet")
    print(chemicals)

    # Predict properties and process the predictions
    cachedir = pathlib.Path('cache') / 'function_cache' / 'chemprop_predictions_cache'
    cachedir.mkdir(parents=True, exist_ok=True)
    pred_func = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)
    for chemical in chemicals:
        llm_input = generate_llm_input(chemical, pred_func)
        

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
