#!/usr/bin/env python3
"""
Script to import a WikiPathways pathway into Cytoscape, attach a numerical property
to each node, color nodes by that property, and export the network image.

Requirements:
  - Cytoscape running with the WikiPathways App installed
  - py4cytoscape, pandas
"""

import argparse
import pandas as pd
import py4cytoscape as p4c

def main():
    parser = argparse.ArgumentParser(
        description="Import a WikiPathways pathway, map a property, and export image."
    )
    parser.add_argument(
        "--pathway", "-p", required=True,
        help="WikiPathways ID (e.g., WP3657) of the pathway to import"
    )
    parser.add_argument(
        "--data", "-d", required=True,
        help="CSV file with columns 'gene' and the specified property"
    )
    parser.add_argument(
        "--property", "-r", dest="prop", required=True,
        help="Name of the numerical property column in the data file"
    )
    parser.add_argument(
        "--vmin", type=float, default=None,
        help="Minimum value for the color gradient (defaults to data min)"
    )
    parser.add_argument(
        "--vmax", type=float, default=None,
        help="Maximum value for the color gradient (defaults to data max)"
    )
    parser.add_argument(
        "--view", choices=["pathway", "network"], default="pathway",
        help="Import as 'pathway' (diagram) or 'network' (topology)"
    )
    args = parser.parse_args()

    # 1) Connect to Cytoscape
    p4c.cytoscape_ping()  # Verify CyREST connectivity :contentReference[oaicite:5]{index=5}

    # 2) Import the pathway from WikiPathways
    cmd = f"wikipathways import-as-{args.view} id={args.pathway}"
    p4c.commands.commands_run(cmd)  # Issue the import command :contentReference[oaicite:6]{index=6}

    # 3) Load the user data
    df = pd.read_csv(args.data)
    df = df[["gene", args.prop]].rename(columns={args.prop: "property"})

    # 4) Merge the data into the Cytoscape node table
    p4c.load_table_data(
        df,
        data_key_column="gene",
        table="node",
        table_key_column="gene"
    )  # Push DataFrame to Cytoscape :contentReference[oaicite:7]{index=7}

    # 5) Create a visual style mapping 'property' → NODE_FILL_COLOR
    vmin = args.vmin if args.vmin is not None else df["property"].min()
    vmax = args.vmax if args.vmax is not None else df["property"].max()
    midpoint = (vmin + vmax) / 2.0

    mapping = p4c.map_visual_property(
        visual_prop="NODE_FILL_COLOR",
        table_column="property",
        mapping_type="continuous",
        table_column_values=[vmin, midpoint, vmax],
        visual_prop_values=["#313695", "#ffffbf", "#a50026"]
    )  # Define color gradient :contentReference[oaicite:8]{index=8}

    style_name = f"Style_{args.prop}"
    p4c.create_visual_style(
        style_name,
        defaults={"NODE_SHAPE": "ELLIPSE", "NODE_SIZE": 30},
        mappings=[mapping]
    )  # Create the style :contentReference[oaicite:9]{index=9}
    p4c.set_visual_style(style_name)  # Apply it :contentReference[oaicite:10]{index=10}

    # 6) Optional layout for network view
    if args.view == "network":
        p4c.layout_network("force-directed")  # Apply force-directed layout :contentReference[oaicite:11]{index=11}

    # 7) Export the final image
    output_file = f"{args.pathway}_{args.prop}.png"
    p4c.export_image(
        filename=output_file,
        type="PNG",
        zoom=200,
        all_graphics_details=True
    )  # Export high‑res image :contentReference[oaicite:12]{index=12}

    print(f"Saved network image to {output_file}")

if __name__ == "__main__":
    main()
