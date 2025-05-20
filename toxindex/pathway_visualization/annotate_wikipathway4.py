#!/usr/bin/env python3
"""
Import a WikiPathways pathway into Cytoscape, map a numeric property onto nodes,
and export a PNG. Supports CSV/TSV and Parquet data inputs, and ensures nodes
without data receive a NaN property (which Cytoscape treats as default).

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
        '--property', '-r', dest='prop', required=True,
        help='Column name in your data file (will be imported as "property")'
    )
    parser.add_argument(
        '--data', '-d', required=True,
        help='CSV, TSV, or Parquet file with at least columns "gene" and your property'
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

    # Connect to Cytoscape
    p4c.cytoscape_ping()

    # Import pathway view
    p4c.commands.commands_run(
        f'wikipathways import-as-{args.view} id={args.pathway}'
    )

    # Load user data
    ext = os.path.splitext(args.data)[1].lower()
    if ext == '.parquet':
        df = pd.read_parquet(args.data)
    else:
        try:
            df = pd.read_csv(args.data, sep=None, engine='python')
        except UnicodeDecodeError:
            df = pd.read_csv(args.data, sep=None, engine='python', encoding='latin1')

    # Validate
    if 'gene' not in df.columns or args.prop not in df.columns:
        sys.exit('ERROR: data must contain columns "gene" and the specified property.')

    # Prepare mapping DataFrame
    df2 = df[['gene', args.prop]].rename(columns={'gene':'name', args.prop:'property'})

    # Load table data into Cytoscape
    p4c.load_table_data(
        df2,
        data_key_column='name',
        table='node',
        table_key_column='name'
    )

    # Define color scale bounds
    vmin = args.vmin if args.vmin is not None else df2['property'].min()
    vmax = args.vmax if args.vmax is not None else df2['property'].max()
    midpoint = (vmin+vmax)/2

    # Create and apply visual style
    style_name = f'Style_{args.prop}'
    defaults = {'NODE_SHAPE':'ELLIPSE','NODE_SIZE':30}
    mapping = p4c.map_visual_property(
        visual_prop='NODE_FILL_COLOR',
        table_column='property',
        mapping_type='continuous',
        table_column_values=[vmin, midpoint, vmax],
        visual_prop_values=['#313695','#ffffbf','#a50026']
    )
    p4c.create_visual_style(style_name, defaults=defaults, mappings=[mapping])
    p4c.set_visual_style(style_name)

    # Optional layout
    if args.view == 'network':
        p4c.layout_network('force-directed')

    # Export image
    out = f"{args.pathway}_{args.prop}.png"
    p4c.export_image(out, type='PNG', resolution=300)
    print(f'Saved image to {out}')

if __name__=='__main__':
    main()
