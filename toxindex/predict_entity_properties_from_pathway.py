from tqdm import tqdm
import pandas as pd
import biobricks as bb
from rdflib import Graph
from rdflib_hdt import HDTStore
import google
from google import genai
from tenacity import retry, wait_fixed, retry_if_exception_type
import argparse
import numpy as np
import os
from dotenv import load_dotenv
import json
import re
from pydantic import BaseModel
from typing import Literal
import pathlib


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


class PathwayEntityPropSchema(BaseModel):
    """Pydantic schema for validating a list of property names returned for a pathway entity."""

    property_names: list[str]


# ──────────────────────────────────────────────────────────────────────────────
# SPARQL helper ─ generic for genes *or* key events
# ──────────────────────────────────────────────────────────────────────────────
def get_pathway_entities(
    identifier: str,
    store: HDTStore,
    kind: Literal["gene", "ke"] = "gene",
) -> pd.DataFrame:
    """
    Return labels for pathway entities:
      • kind="gene" → WikiPathways GeneProducts
      • kind="ke"   → AOP‑Wiki Key Events
    """
    g = Graph(store=store)

    if kind == "gene":                    # ── WikiPathways branch ─────────────
        q = f"""
        PREFIX wp:   <http://vocabularies.wikipathways.org/wp#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX dct:  <http://purl.org/dc/terms/>

        SELECT DISTINCT ?gp ?label WHERE {{
            ?pw  a wp:Pathway ;
                 dct:identifier "{identifier}" .
            ?gp  a wp:GeneProduct ;
                 dct:isPartOf ?pw ;
                 rdfs:label ?label .
        }}
        ORDER BY ?label
        """
        cols = ["geneProduct", "geneLabel"]

    elif kind == "ke":                    # ── AOP‑Wiki branch ─────────────────
        # If the caller already passed the full IRI, keep it; otherwise build it
        if identifier.startswith("http"):
            aop_uri = identifier
        else:
            aop_uri = f"https://identifiers.org/aop/{identifier}"

        q = f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX dc:   <http://purl.org/dc/elements/1.1/>
        PREFIX aopo: <http://aopkb.org/aop_ontology#>
        PREFIX obo:  <http://purl.obolibrary.org/obo/>

        SELECT DISTINCT ?ke (COALESCE(?l1,?l2) AS ?label) WHERE {{
            VALUES ?aop {{ <{aop_uri}> }}
            VALUES ?pred {{ aopo:has_key_event aopo:hasKeyEvent obo:aopo_0001015 }}

            ?aop ?pred ?ke .

            # accept both AOPO namespace variants for the KE class
            VALUES ?keClass {{ aopo:KeyEvent  obo:aopo_0000044 }}
            ?ke a ?keClass .

            OPTIONAL {{ ?ke rdfs:label ?l1 }}
            OPTIONAL {{ ?ke dc:title   ?l2 }}
        }}
        ORDER BY ?label
        """
        cols = ["keyEvent", "eventLabel"]
    
    return pd.DataFrame(g.query(q), columns=cols)


# def get_pathway_entities(identifier: str, store: HDTStore,
#                          kind: Literal["gene", "ke"] = "gene") -> pd.DataFrame:
#     """Return labels for pathway entities (gene products or AOP key events)."""

#     if kind == "gene":                           # ── WikiPathways branch ──
#         q = f"""
#         PREFIX wp:   <http://vocabularies.wikipathways.org/wp#>
#         PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
#         PREFIX dct:  <http://purl.org/dc/terms/>

#         SELECT DISTINCT ?gp ?label WHERE {{
#             ?pw  a wp:Pathway ;
#                  dct:identifier "{identifier}" .
#             ?gp  a wp:GeneProduct ;
#                  dct:isPartOf ?pw ;
#                  rdfs:label ?label .
#         }}
#         ORDER BY ?label
#         """
#         cols = ["geneProduct", "geneLabel"]

#     else:                                        # ── AOP‑Wiki branch ──────
#         q = f"""
#         PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
#         PREFIX dc:   <http://purl.org/dc/elements/1.1/>
#         PREFIX dct:  <http://purl.org/dc/terms/>
#         PREFIX aopo: <http://aopkb.org/aop_ontology#>
#         PREFIX obo:  <http://purl.obolibrary.org/obo/>

#         SELECT DISTINCT ?ke (COALESCE(?l1,?l2) AS ?label) WHERE {{
#             # -------- pathway node -------------------------------------------------
#             VALUES ?aopClass {{ aopo:AdverseOutcomePathway  obo:aopo_0000001 }}
#             ?aop  a ?aopClass ;
#                   (dc:identifier|dct:identifier) "{identifier}" ;
#                   (aopo:has_key_event|aopo:hasKeyEvent|obo:aopo_0001015) ?ke .

#             # -------- key‑event node -----------------------------------------------
#             VALUES ?keClass {{ aopo:KeyEvent  obo:aopo_0000044 }}
#             ?ke   a ?keClass .
#             OPTIONAL {{ ?ke rdfs:label ?l1 }}
#             OPTIONAL {{ ?ke dc:title   ?l2 }}
#         }}
#         ORDER BY ?label
#         """
#         cols = ["keyEvent", "eventLabel"]

#     g = Graph(store=store)
#     return pd.DataFrame(g.query(q), columns=cols)


# ──────────────────────────────────────────────────────────────────────────────
# Utility for sloppy JSON returned by LLM
# ──────────────────────────────────────────────────────────────────────────────

def get_response_json(response):
    """Attempt to coerce an LLM response into valid JSON."""

    raw = response if isinstance(response, str) else response.text

    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        snippet = re.search(r"\{.*", raw, re.S).group(0)
        if snippet.count("\"") % 2 == 1:
            snippet += "\""
        if snippet.count("[") > snippet.count("]"):
            snippet += "]"
        if snippet.count("{") > snippet.count("}"):
            snippet += "}"
        return json.loads(snippet)


# ──────────────────────────────────────────────────────────────────────────────
# Gemini wrapper
# ──────────────────────────────────────────────────────────────────────────────

@retry(
    retry=retry_if_exception_type(google.genai.errors.ServerError),
    wait=wait_fixed(10),
    reraise=True,
)
def get_property_response(client, entity_label: str, pathway: str, properties: list[str]):
    """Call Gemini with automatic retries on server‐side errors."""

    prompt = (
        "Select properties from the list below that are strongly associated "
        f"with the pathway entity {entity_label} in the biological pathway {pathway}.\n\n"
        f"{properties}\n\n"
        "Output one unique property per line and nothing else; copy the property names exactly."
    )

    return client.models.generate_content(
        model="gemini-2.0-flash",
        contents=[prompt],
        config={"temperature": 0.0},
    )


# ──────────────────────────────────────────────────────────────────────────────
# LLM‑driven property selection and downstream processing
# ──────────────────────────────────────────────────────────────────────────────

def predict_properties_for_entities(
    outdir: pathlib.Path,
    client,
    pathway: str,
    properties: list[str],
    pathway_entities: pd.DataFrame,
    prefix: str,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Ask Gemini which *properties* apply to each pathway entity."""

    fname = outdir / f"{pathway}.{prefix}_property_predictions.parquet"
    if fname.exists() and use_cache:
        print(f"Using cached file: {fname}")
        return pd.read_parquet(fname)

    rows = []
    for entity_label in tqdm(
        pathway_entities.iloc[:, 1], desc="Predicting properties", unit="entity"
    ):
        try:
            response = get_property_response(client, entity_label, pathway, properties)
        except Exception as exc:
            print(f"Failed for entity {entity_label!r}: {exc}")
            continue

        prop_list = [p for p in response.text.splitlines() if p in properties]
        rows.append({"entity": entity_label, "pathway": pathway, "property_list": prop_list})

    df = pd.DataFrame(rows)
    df.to_parquet(fname)
    return df


def entity_chemical_strengths(
    outdir: pathlib.Path,
    pathway: str,
    entity_property_predictions: pd.DataFrame,
    predictions: pd.DataFrame,
    prefix: str,
    use_cache: bool = True,
) -> pd.DataFrame:
    """Filter chemical predictions by the property list for each entity."""

    fname = outdir / f"{pathway}.{prefix}_property_chemicals.parquet"
    if fname.exists() and use_cache:
        df = pd.read_parquet(fname)
        df["chemicals_values"] = df["chemicals_values"].apply(json.loads)
        print(f"Using cached file: {fname}")
        return df

    rows = []
    for _, epp in tqdm(
        entity_property_predictions.iterrows(),
        desc="Selecting chemicals",
        total=len(entity_property_predictions),
    ):
        chems = [
            (pred["inchi"], pred["value"])
            for _, pred in predictions.iterrows()
            if pred["property_title"] in epp["property_list"]
        ]
        rows.append({"entity": epp["entity"], "pathway": pathway, "chemicals_values": chems})

    df = pd.DataFrame(rows)
    df["chemicals_values"] = df["chemicals_values"].apply(json.dumps)
    df.to_parquet(fname)
    df["chemicals_values"] = df["chemicals_values"].apply(json.loads)
    return df


def softmax(arr, beta=1.0):
    """Compute the softmax of an array."""
    e_x = np.exp(beta*arr)
    weights = e_x / np.sum(e_x)
    return np.dot(arr, weights)


# ──────────────────────────────────────────────────────────────────────────────
# Main driver
# ──────────────────────────────────────────────────────────────────────────────

def main(
    cachedir: pathlib.Path,
    projectdir: pathlib.Path,
    outdir: pathlib.Path,
    pathway: str,
    kind: Literal["gene", "ke"],
    use_cache_predictions: bool = True,
    use_cache_chemicals: bool = True,
    # softmax_beta: float = 1.0,
):
    load_dotenv()
    client = genai.Client(api_key=os.environ["GEMINI_API_KEY"])

    predictions = pd.read_parquet(projectdir / "predict_chemicals" / "chemprop_predictions.parquet")

    if kind == "gene":
        wikipathways = bb.assets("wikipathways")
        store = HDTStore(wikipathways.wikipathways_hdt)
    elif kind == "ke":
        aopwiki = bb.assets("aopwikirdf-kg")
        store = HDTStore(aopwiki.AOPWikiRDF_hdt)

    print("Retrieving pathway entities...")
    pathway_entities = get_pathway_entities(pathway, store, kind=kind)
    store.close()

    if (pathway_entities is None) or pathway_entities.empty:
        raise RuntimeError("No entities found for the specified pathway.")

    properties = [p.strip() for p in (cachedir / "predicted_property_names.txt").read_text().splitlines()]

    prefix = kind  # "gene" or "ke" – keeps filenames readable

    entity_property_predictions = predict_properties_for_entities(
        outdir,
        client,
        pathway,
        properties,
        pathway_entities,
        prefix,
        use_cache_predictions,
    )

    entity_property_chemicals = entity_chemical_strengths(
        outdir,
        pathway,
        entity_property_predictions,
        predictions,
        prefix,
        use_cache_chemicals,
    )

    merged = entity_property_chemicals.groupby("entity", as_index=False).agg({
        "chemicals_values": lambda lsts: sum(lsts, []),
    })

    def compute_stats(chem_vals):
        values = [v for _, v in chem_vals]
        if not values:
            return pd.Series(dtype=float)
        arr = np.array(values)
        return pd.Series(
            {
                "sum_value": arr.sum(),
                "mean_value": arr.mean(),
                # "2-norm": np.linalg.norm(arr),
                # "3-norm": np.linalg.norm(arr, ord=3),
                # "4-norm": np.linalg.norm(arr, ord=4),
                "max_value": arr.max(),
                "softmax" : softmax(arr, beta = 1.0),
                "softmax2": softmax(arr, beta = 2.0),
                "softmax3": softmax(arr, beta = 3.0),
            }
        )

    summary = pd.concat([merged.drop(columns="chemicals_values"), merged["chemicals_values"].apply(compute_stats)], axis=1)

    fname_summary = outdir / f"{pathway}.{prefix}_property_chemicals_summary.parquet"
    summary.to_parquet(fname_summary)
    print(summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Property extraction for genes or key events in a pathway.")
    parser.add_argument("--project", type=str, default="hepatotoxic", help="Project name")
    parser.add_argument("--pathway", type=str, required=True, help="Pathway ID to query")
    parser.add_argument("--kind", choices=["gene", "ke"], default="gene", help="Entity type: gene (WikiPathways) or ke (AOP‑Wiki)")
    parser.add_argument("--use_cache_predictions", type=str2bool, default=True, help="Reuse cached entity→property predictions")
    parser.add_argument("--use_cache_chemicals", type=str2bool, default=True, help="Reuse cached chemicals filtered by entity properties")
    # parser.add_argument("--softmax_beta", type=float, default=1.0, help="Softmax β for weighted averages")
    args = parser.parse_args()

    cachedir = pathlib.Path("cache")
    projectdir = cachedir / "projects" / args.project
    outdir = projectdir / f"{args.kind}_property_predictions"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"cachedir = {cachedir}, projectdir = {projectdir}, outdir = {outdir}")

    main(
        cachedir,
        projectdir,
        outdir,
        args.pathway,
        args.kind,
        args.use_cache_predictions,
        args.use_cache_chemicals,
        # args.softmax_beta,
    )
