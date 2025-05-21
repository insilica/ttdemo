#!/usr/bin/env python3
"""
Import a WikiPathways pathway (or an AOP-Wiki network) into Cytoscape, map a
numeric property onto its nodes, and export a PNG.  The script now works for
**genes** (WikiPathways gene products) *or* **key events** (AOP-Wiki), selected
with ``--kind gene|ke``.

Changes from the gene-only version are limited to:
  • replacing hard-coded references to "gene" with the neutral term *entity* in
    data handling, while
  • keeping the ``gene`` / ``ke`` token in filenames so downstream tooling is
    unaffected.

Requirements (unchanged):
  - Cytoscape (>=3.8) running with WikiPathways app + cyREST
  - py4cytoscape  •  pandas  •  requests (indirect)
"""

import argparse
import pathlib
import sys

import pandas as pd
import py4cytoscape as p4c

import tempfile
import io
import os
import requests

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


def xml_to_graphml(xml_bytes: bytes) -> bytes:
    root = ET.fromstring(xml_bytes)
    G = nx.DiGraph()

    # 1) Extract Key Events
    for ke in root.findall('.//{*}key-event'):
        ke_id    = ke.get('id')
        if not ke_id:
            continue
        ke_title = ke.findtext('{*}title', default='')
        G.add_node(ke_id, label=ke_title)

    # 2) Extract Key Event Relationships
    for ker in root.findall('.//{*}key-event-relationship'):
        src = ker.findtext('{*}ke-upstream-id')
        tgt = ker.findtext('{*}ke-downstream-id')
        rel = ker.findtext('{*}relationship-type', default='')
        if src and tgt:
            G.add_edge(src, tgt, rel=rel)

    # 3) (Optional) Debug counts
    print(f"Built graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

    # 4) Write GraphML
    buf = io.BytesIO()
    nx.write_graphml(G, buf, encoding='utf-8')
    return buf.getvalue()


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
            xml = fetch_aop_xml(int(args.pathway))
            gml_bytes = xml_to_graphml(xml)
            import_graphml_bytes(gml_bytes)
            # import_graphml_bytes(gml_bytes, collection="AOP-Wiki")
            # p4c.import_network_from_string(gml.decode("utf-8"), file_type="graphml")
        except Exception as e:
            raise RuntimeError(f"could not import AOP {args.pathway}: {e}")
    else:
        raise ValueError(f"Unknown entity type: {args.kind}")

    # ── 3) Locate data file ─────────────────────────────────────────────────
    project_dir = pathlib.Path("cache") / "projects" / args.project
    if args.data is None:
        args.data = (
            project_dir
            / f"{args.kind}_property_predictions"
            / f"{args.pathway}.{args.kind}_property_chemicals_summary.parquet"
        )
    args.data = pathlib.Path(args.data)

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
    df2 = df.loc[:, ["entity", args.prop]].copy().rename(columns={"entity": "name", args.prop: "property"})

    # ── 7) Push to Cytoscape node table ──────────────────────────────────────
    p4c.load_table_data(df2, data_key_column="name", table="node", table_key_column="name")

    # ── 8) Determine colour-scale bounds ─────────────────────────────────────
    vmin = args.vmin if args.vmin is not None else df2["property"].min()
    vmax = args.vmax if args.vmax is not None else df2["property"].max()
    midpoint = (vmin + vmax) / 2.0

    # ── 9) Create visual style ───────────────────────────────────────────────
    style_name = f"Style_{args.prop}"
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
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
        table_column="name",
        mapping_type="passthrough",
    )
    p4c.create_visual_style(style_name, defaults=defaults, mappings=[fill_mapping, label_mapping])
    p4c.set_visual_style(style_name)

    # ── 10) Optional layout for topology view ────────────────────────────────
    if args.view == "network":
        p4c.layout_network("force-directed")

    # ── 11) Export image ─────────────────────────────────────────────────────
    image_dir = project_dir / "images"
    image_dir.mkdir(parents=True, exist_ok=True)
    out_png = image_dir / f"{args.pathway}_{args.prop}.png"
    p4c.export_image(str(out_png), type="PNG", resolution=300)
    print(f"Saved image to {out_png}")


if __name__ == "__main__":
    main()
