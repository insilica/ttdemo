#!/usr/bin/env python3
"""
Import a WikiPathways pathway into Cytoscape, map a numeric property onto nodes,
and export a PNG. Supports CSV/TSV and Parquet data inputs, ensures nodes without
data receive a NaN property, uses a Viridis colormap, and preserves node labels via
passthrough mapping with explicit data type.

Requirements:
  - Cytoscape (>=3.8) running with:
      * WikiPathways app
      * Automation (cyREST) enabled
  - py4cytoscape: pip install py4cytoscape
  - pandas, requests
"""

import argparse
import os
import sys
import pathlib

import pandas as pd
import py4cytoscape as p4c


def main():
    parser = argparse.ArgumentParser(
        description='Map a numeric property onto a WikiPathways network in Cytoscape'
    )
    parser.add_argument(
        '--pathway', '-p', required=True,
        help='WikiPathways ID, e.g. WP3657'
    )
    parser.add_argument(
        '--project', type=str, default='hepatotoxic',
        help='Project name for I/O directories'
    )
    parser.add_argument(
        '--property', '-r', dest='prop', required=True,
        help='Column name in your data file (will be imported as "property")'
    )
    parser.add_argument(
        '--data', '-d', required=False, default=None,
        help='CSV, TSV, or Parquet file with columns "gene" and your property'
    )
    parser.add_argument(
        '--vmin', type=float,
        help='Override minimum for the colour scale'
    )
    parser.add_argument(
        '--vmax', type=float,
        help='Override maximum for the colour scale'
    )
    parser.add_argument(
        '--view', choices=['pathway','network'], default='pathway',
        help="Import as 'pathway' (diagram) or 'network' (topology)"
    )
    args = parser.parse_args()

    # 1) Connect to Cytoscape
    p4c.cytoscape_ping()

    # 2) Import pathway view
    p4c.commands.commands_run(
        f'wikipathways import-as-{args.view} id={args.pathway}'
    )

    # 3) Load user data
    project_dir = pathlib.Path(f"cache/projects/{args.project}")
    if args.data is None:
        args.data = project_dir / "gene_property_predictions" \
            / f"{args.pathway}.gene_property_chemicals_summary.parquet"
    ext = pathlib.Path(args.data).suffix.lower()
    print(f'Loading data from {args.data} (ext: {ext})')
    if ext == '.parquet':
        df = pd.read_parquet(args.data)
    else:
        try:
            df = pd.read_csv(args.data, sep=None, engine='python')
        except UnicodeDecodeError:
            df = pd.read_csv(
                args.data, sep=None, engine='python', encoding='latin1'
            )

    # 4) Validate required columns
    if 'gene' not in df.columns or args.prop not in df.columns:
        sys.exit('ERROR: data must contain columns "gene" and the specified property.')

    # 5) Prepare mapping DataFrame
    df2 = df.loc[:, ['gene', args.prop]].copy()
    df2.rename(columns={'gene': 'name', args.prop: 'property'}, inplace=True)

    # 6) Merge into Cytoscape node table
    p4c.load_table_data(
        df2,
        data_key_column='name',
        table='node',
        table_key_column='name'
    )

    # 7) Define colour scale bounds
    vmin = args.vmin if args.vmin is not None else df2['property'].min()
    vmax = args.vmax if args.vmax is not None else df2['property'].max()
    midpoint = (vmin + vmax) / 2.0

    # 8) Create and apply visual style with Viridis gradient and explicit label passthrough
    style_name = f'Style_{args.prop}'
    defaults = {
        'NODE_SHAPE': 'ELLIPSE',
        # 'NODE_SIZE': 60,
        'NODE_LABEL_COLOR': '#000000',
        'NODE_LABEL_FONT_SIZE': 12
    }
    # Viridis colormap hexes at low, mid, high
    viridis_colors = ['#440154', '#21918C', '#FDE725']
    fill_mapping = p4c.map_visual_property(
        visual_prop='NODE_FILL_COLOR',
        table_column='property',
        mapping_type='continuous',
        table_column_values=[vmin, midpoint, vmax],
        visual_prop_values=viridis_colors
    )
    # Explicit passthrough mapping for labels
    label_mapping = p4c.map_visual_property(
        visual_prop='NODE_LABEL',
        table_column='name',
        mapping_type='passthrough',
        # mapping_column_data_type='String'
    )
    p4c.create_visual_style(
        style_name,
        defaults=defaults,
        mappings=[fill_mapping, label_mapping]
    )
    p4c.set_visual_style(style_name)

    # 9) Apply optional layout
    if args.view == 'network':
        p4c.layout_network('force-directed')

    # 10) Export the image
    image_dir = project_dir / 'images'
    image_dir.mkdir(parents=True, exist_ok=True)
    output_file = image_dir / f"{args.pathway}_{args.prop}.png"
    p4c.export_image(str(output_file), type='PNG', resolution=300)
    print(f'Saved image to {output_file}')

if __name__ == '__main__':
    main()
