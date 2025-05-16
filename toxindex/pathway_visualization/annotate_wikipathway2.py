#!/usr/bin/env python3
"""
Import a WikiPathways pathway into Cytoscape, map a numeric property to node colours,
and export a PNG.

Requirements:
  - Cytoscape (>=3.8) running with:
      * WikiPathways app
      * Automation (cyREST) enabled
  - py4cytoscape Python package: pip install py4cytoscape
"""

import argparse
import os
import sys

import pandas as pd
import requests
import py4cytoscape as p4c


def main():
    p = argparse.ArgumentParser(
        description='Map a numeric property onto a WikiPathways network in Cytoscape'
    )
    p.add_argument(
        '--pathway', '-p',
        required=True,
        help='WikiPathways ID, e.g. WP3657'
    )
    p.add_argument(
        '--property', '-r', dest='prop',
        required=True,
        help='Column name in your data file (will be imported as "property")'
    )
    p.add_argument(
        '--data', '-d',
        required=True,
        help='Parquet file with at least columns "gene" and your property'
    )
    p.add_argument(
        '--vmin', type=float,
        help='Override minimum for the colour scale'
    )
    p.add_argument(
        '--vmax', type=float,
        help='Override maximum for the colour scale'
    )
    args = p.parse_args()

    # 1) Fetch GPML
    gpml_url = (
        f'https://webservice.wikipathways.org/getPathway'
        f'?pathwayId={args.pathway}&format=gpml'
    )
    print(f'Downloading GPML from {gpml_url}...')
    # r = requests.get(gpml_url)
    # r.raise_for_status()
    # gpml_file = f'{args.pathway}.gpml'
    # with open(gpml_file, 'wb') as fh:
    #     fh.write(r.content)
    # print(f'  → saved to {gpml_file}')
    gpml_file = f'toxindex/pathway_visualization/{args.pathway}.gpml'

    # 2) Import into Cytoscape
    print('Connecting to Cytoscape...')
    p4c.cytoscape_ping()
    print('Importing network into Cytoscape...')
    p4c.import_network_from_file(gpml_file)

    # 3) Read your property table and import as node table
    print(f'Loading data from {args.data}...')
    df = pd.read_parquet(args.data)
    if 'gene' not in df.columns or args.prop not in df.columns:
        raise ValueError('Data file must contain "gene" and the specified property column.')

    df2 = df[['gene', args.prop]].rename(columns={args.prop: 'property'})
    tmp_csv = f'tmp_{args.pathway}_{args.prop}.csv'
    df2.to_csv(tmp_csv, index=False)
    print('Importing property values into Cytoscape node table...')
    p4c.commands.commands_call(
        f'table import file file="{tmp_csv}" table="node" keyColumnLabel="gene"'
    )
    os.remove(tmp_csv)

    # 4) Build & apply a continuous colour style
    vmin = args.vmin if args.vmin is not None else df2['property'].min()
    vmax = args.vmax if args.vmax is not None else df2['property'].max()
    style_name = f'Style_{args.prop}'
    defaults = {
        'NODE_SHAPE': 'ELLIPSE',
        'NODE_SIZE': 30
    }
    mapping = {
        'mappingType': 'continuous',
        'mappingColumn': 'property',
        'mappingColumnDataType': 'Double',
        'visualProperty': 'NODE_FILL_COLOR',
        'points': [
            { 'value': vmin, 'lesser': '#313695', 'equal': '#313695', 'greater': '#313695' },
            { 'value': (vmin + vmax)/2, 'lesser': '#ffffbf', 'equal': '#ffffbf', 'greater': '#ffffbf' },
            { 'value': vmax, 'lesser': '#a50026', 'equal': '#a50026', 'greater': '#a50026' }
        ]
    }

    print(f'Creating visual style "{style_name}" (range {vmin}–{vmax})...')
    p4c.create_visual_style(style_name, defaults=defaults, mappings=[mapping])
    p4c.set_visual_style(style_name)

    # 5) Layout and finalise
    print('Applying force-directed layout...')
    p4c.layout_network('force-directed')

    # 6) Export
    out_png = f"{args.pathway}_{args.prop}.png"
    print(f'Exporting to {out_png}...')
    p4c.export_image(out_png, type='PNG', resolution=300)
    print('Done.')


if __name__ == '__main__':
    main()
