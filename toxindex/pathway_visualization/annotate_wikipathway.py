#!/usr/bin/env python3
"""
Import a WikiPathways pathway (or an AOP-Wiki network) into Cytoscape, map a
Import a WikiPathways pathway (or an AOP-Wiki network) into Cytoscape, map a
numeric property onto its nodes, and export a PNG.  The script now works for
**genes** (WikiPathways gene products) *or* **key events** (AOP-Wiki), selected
**genes** (WikiPathways gene products) *or* **key events** (AOP-Wiki), selected
with ``--kind gene|ke``.

Changes from the gene-only version are limited to:
  • replacing hard-coded references to "gene" with the neutral term *entity* in
Changes from the gene-only version are limited to:
  • replacing hard-coded references to "gene" with the neutral term *entity* in
    data handling, while
  • keeping the ``gene`` / ``ke`` token in filenames so downstream tooling is
    unaffected.

Requirements (unchanged):
  - Cytoscape (>=3.8) running with WikiPathways app + cyREST
  - py4cytoscape  •  pandas  •  requests (indirect)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
import py4cytoscape as p4c

import tempfile
import os
import requests
from typing import Union, Dict, Set, List

import xml.etree.ElementTree as ET
import networkx as nx
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def fetch_aop_xml(aop_id: int) -> bytes:
    """
    Download the AOP-Wiki XML feed, ignoring the server’s incomplete
    certificate chain.  If that still fails, try the plain-HTTP mirror.
    """
    urls = [
        f"https://aopwiki.org/aops/{aop_id}.xml",   # official (bad TLS chain)
        f"http://aopwiki.org/aops/{aop_id}.xml",    # falls back to 302→https
        f"https://aopwiki-api.bio.tools/aops/{aop_id}.xml",  # mirror, good TLS
    ]
    for u in urls:
        try:
            r = requests.get(u, timeout=20, verify=False)  # <-- key change
            r.raise_for_status()
            if r.text.lstrip().startswith("<?xml"):
                return r.content
        except Exception as e:
            print(f"[WARN] {u} → {e.__class__.__name__}: {e}")
    raise RuntimeError(f"Unable to fetch XML for AOP {aop_id}")


def aop_xml_bytes_to_graphml(xml_bytes: bytes) -> Path:
    """
    Convert AOP-Wiki XML (as bytes) into a temporary GraphML file and
    return the path to that file.
    """

    # --------------------------- parse XML ----------------------------------
    root = ET.fromstring(xml_bytes)
    ns_url = root.tag.partition('}')[0].lstrip('{')
    ns = {"a": ns_url}  # single dynamic namespace mapping

    # --------------------- gather global lookup tables ----------------------
    aop_ids: List[str] = sorted({e.get("aop-wiki-id")
                                 for e in root.findall(".//a:aop-reference", ns)
                                 if e.get("aop-wiki-id")}) or [""]

    ke_wiki_id: Dict[str, str] = {
        ref.get("id"): ref.get("aop-wiki-id", "?")
        for ref in root.findall(".//a:key-event-reference", ns)
    }

    mie_ids: Set[str] = {
        e.get("key-event-id") or e.get("key_event_id")
        for e in root.findall(".//a:molecular-initiating-event", ns)
    }
    ao_ids: Set[str] = {
        e.get("key-event-id") or e.get("key_event_id")
        for e in root.findall(".//a:adverse-outcome", ns)
    }

    # --------------------------- build graph --------------------------------
    G = nx.DiGraph()

    # ---- nodes
    for ke in root.findall(".//a:key-event", ns):
        ke_id = ke.get("id")
        if not ke_id:
            continue
        title = ke.findtext("a:title", default="", namespaces=ns)

        if ke_id in mie_ids:
            ke_type = "Molecular Initiating Event"
        elif ke_id in ao_ids:
            ke_type = "Adverse Outcome"
        else:
            ke_type = "Key Event"

        G.add_node(
            ke_id,
            **{
                "ke_in_aop": str(aop_ids),
                "ke_type": ke_type,
                "label": f"KE {ke_wiki_id.get(ke_id, ke_id)}",
                "name": title,
                "selected": "false",
                "shared name": title,
                "value": title,
            },
        )

    # ---- edges (handle several schema spellings)
    for ker in root.findall(".//a:key-event-relationship", ns):
        src = (
            ker.findtext(".//a:upstream-id", namespaces=ns)
            or ker.findtext(".//a:upstream_id", namespaces=ns)
            or ker.findtext(".//a:ke-upstream-id", namespaces=ns)
            or ker.findtext(".//a:ke_upstream_id", namespaces=ns)
        )
        dst = (
            ker.findtext(".//a:downstream-id", namespaces=ns)
            or ker.findtext(".//a:downstream_id", namespaces=ns)
            or ker.findtext(".//a:ke-downstream-id", namespaces=ns)
            or ker.findtext(".//a:ke_downstream_id", namespaces=ns)
        )
        if src and dst and src in G and dst in G:
            G.add_edge(src, dst)

    # ------------------------- write GraphML --------------------------------
    with tempfile.NamedTemporaryFile(
        mode="w+", delete=False, suffix=".graphml", encoding="utf-8"
    ) as fh:
        graphml_path = Path(fh.name)
        nx.write_graphml(G, graphml_path)

    print(
        f"[INFO] Wrote {len(G.nodes)} nodes and {len(G.edges)} edges → {graphml_path}"
    )
    return graphml_path


def import_graphml_bytes(gml_bytes: bytes, *, collection=None):
    """Load raw GraphML bytes into Cytoscape (all py4cytoscape versions)."""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".graphml") as fh:
        fh.write(gml_bytes)
        tmp_path = fh.name
    try:
        p4c.import_network_from_file(tmp_path)   # ← works in every build :contentReference[oaicite:0]{index=0}
        # p4c.import_network_from_file(tmp_path, collection=collection)   # ← works in every build :contentReference[oaicite:0]{index=0}
    finally:
        os.unlink(tmp_path)


def load_aop_graphml(aop_id: int, *, collection=None):
    url = ("https://aopwiki-api.cloud.vhp4safety.nl/"
           "api-git/marvinm2/AOPWikiQueries/get-aop-graphml")
    params  = {"aop_id": aop_id, "layout": "true"}
    r = requests.get(url, params=params,
                     headers={"accept": "application/graphml+xml"},
                     timeout=30)
    r.raise_for_status()                              # 200 → GraphML bytes
    with tempfile.NamedTemporaryFile(delete=False, suffix=".graphml") as fh:
        fh.write(r.content); tmp = fh.name
    try:
        p4c.import_network_from_file(tmp, collection=collection)
    finally:
        os.unlink(tmp)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map a numeric property onto a pathway/network in Cytoscape",
    )
    parser.add_argument("--pathway", "-p", required=True, help="WikiPathways or AOP-Wiki ID, e.g. WP3657 or 37")
    parser.add_argument("--project", default="hepatotoxic", help="Project name for I/O directories")
    parser.add_argument("--property", "-r", dest="prop", required=True, help="Column name to visualise")
    parser.add_argument("--kind", choices=["gene", "ke"], default="gene", help="Entity type contained in the data file")
    parser.add_argument("--data", "-d", default=None, help="CSV/TSV/Parquet with columns 'entity' and property column")
    parser.add_argument("--vmin", type=float, help="Override colour-scale minimum")
    parser.add_argument("--vmax", type=float, help="Override colour-scale maximum")
    parser.add_argument("--view", choices=["pathway", "network"], default="pathway", help="Import as diagram or topology view")
    args = parser.parse_args()

    # ── 1) Ping Cytoscape ────────────────────────────────────────────────────
    p4c.cytoscape_ping()

    # ── 2) Import pathway / network ──────────────────────────────────────────
    if args.kind == "gene":
        # WikiPathways
        p4c.commands.commands_run(f"wikipathways import-as-{args.view} id={args.pathway}")
    elif args.kind == "ke":
        # AOP-Wiki
        try:
            xml_bytes = fetch_aop_xml(int(args.pathway))
            graphml_path = aop_xml_bytes_to_graphml(xml_bytes)
            p4c.import_network_from_file(graphml_path)
            graphml_path.unlink()  # remove temp file
            # gml_bytes = xml_to_graphml(xml)
            # import_graphml_bytes(gml_bytes)

            # # hardcoded workaround for temporary inspection
            # p4c.import_network_from_file('/home/john/Downloads/AOP_21052025.graphml')
        except Exception as e:
            raise RuntimeError(f"could not import AOP {args.pathway}: {e}")
    else:
        raise ValueError(f"Unknown entity type: {args.kind}")

    # ── 3) Locate data file ─────────────────────────────────────────────────
    project_dir = Path("cache") / "projects" / args.project
    if args.data is None:
        args.data = (
            project_dir
            / f"{args.kind}_property_predictions"
            / f"{args.pathway}.{args.kind}_property_chemicals_summary.parquet"
        )
    args.data = Path(args.data)

    # ── 4) Load data ─────────────────────────────────────────────────────────
    print(f"Loading data from {args.data}")
    if args.data.suffix.lower() == ".parquet":
        df = pd.read_parquet(args.data)
    else:
        df = pd.read_csv(args.data, sep=None, engine="python", encoding="utf8", on_bad_lines="skip")

    # ── 5) Validate required columns ─────────────────────────────────────────
    if "entity" not in df.columns or args.prop not in df.columns:
        sys.exit("ERROR: data must contain columns 'entity' and the specified property.")

    # ── 6) Prepare mapping frame ─────────────────────────────────────────────
    col_label = "name" if args.kind == "gene" else "label"
    df2 = df.loc[:, ["entity", args.prop]].copy().rename(columns={"entity": col_label, args.prop: "property"})

    # ── 7) Push to Cytoscape node table ──────────────────────────────────────
    p4c.load_table_data(df2, data_key_column=col_label, table="node", table_key_column=col_label)

    # ── 8) Determine colour-scale bounds ─────────────────────────────────────
    vmin = args.vmin if args.vmin is not None else df2["property"].min()
    vmax = args.vmax if args.vmax is not None else df2["property"].max()
    midpoint = (vmin + vmax) / 2.0

    # ── 9) Create visual style ───────────────────────────────────────────────
    style_name = f"Style_{args.prop}"
    defaults = {
        # "NODE_SHAPE": "ELLIPSE",
        "NODE_SHAPE": "RECTANGLE",
        "NODE_LABEL_COLOR": "#000000",
        "NODE_LABEL_FONT_SIZE": 12,
    }
    viridis = ["#440154", "#21918C", "#FDE725"]  # low → mid → high
    fill_mapping = p4c.map_visual_property(
        visual_prop="NODE_FILL_COLOR",
        table_column="property",
        mapping_type="continuous",
        table_column_values=[vmin, midpoint, vmax],
        visual_prop_values=viridis,
    )
    label_mapping = p4c.map_visual_property(
        visual_prop="NODE_LABEL",
        # table_column="name" if args.kind == "gene" else "label",
        table_column="name",
        mapping_type="passthrough",
    )
    p4c.create_visual_style(style_name, defaults=defaults, mappings=[fill_mapping, label_mapping])
    p4c.set_visual_style(style_name)

    # ── 10) Optional layout for graph view ───────────────────────────────────
    if args.view == "network":
        p4c.layout_network("force-directed")
    else:
        p4c.layout_network("hierarchical")

    # ── 11) Export image ─────────────────────────────────────────────────────
    image_dir = project_dir / "images"
    image_dir.mkdir(parents=True, exist_ok=True)
    out_png = image_dir / f"{args.pathway}_{args.prop}.png"
    p4c.export_image(str(out_png), type="PNG", resolution=300)
    print(f"Saved image to {out_png}")


if __name__ == "__main__":
    main()
