import pubchempy as pcp
import os
from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors
import time
import sys
import warnings
import pathlib
import pandas as pd
from toxindex.utils.helper import handle_exceptions, rate_limit_lockfile
import logging
import toxindex.utils.simplecache as simplecache
from tqdm import tqdm
import requests

# Suppress RDKit warnings for cleaner output (optional)
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')
warnings.filterwarnings("ignore", category=UserWarning, module='rdkit')

def parse_chemicals(input_path, output_path):
    """
    Parse chemical names from a file and retrieve their information from PubChem.
    
    Args:
        input_path (str): Path to the file containing chemical names
        output_path (str): Path to save the output CSV file
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    logging.info(f"Processing chemicals from {input_path}")
    
    # Read chemical names from input file
    chemical_names = set(line.strip() for line in pathlib.Path(input_path).read_text().splitlines() if line.strip())
    chemical_names = sorted(list(chemical_names))
    
    results = []
    sc = simplecache.simple_cache(pathlib.Path('cache/function_cache/parse_chemicals'))

    def try_alternative_name(query,max_results=5):
        esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'pccompound',
            'term': query,
            'retmode': 'json',
            'retmax': max_results
        }
        r = requests.get(esearch_url, params=params)
        cids = r.json().get('esearchresult').get('idlist')

        records = []
        for cid in cids:
            try:
                compound = pcp.get_compounds(cid, 'cid')[0]
                if compound.inchi:
                    records.append((compound, compound.inchi))
            except Exception as e:
                continue

        if not records:
            return None

        # Sort by length of InChI
        records.sort(key=lambda x: len(x[1]))
        best = records[0][0]
        return {"name": query, "cid": best.cid, "inchi": best.inchi}
    
    def parse_chemical(name):
        logging.info(f"Parsing chemical: {name}")
        time.sleep(0.33)
        compounds = pcp.get_compounds(name, 'name')

        # if compound found, return info
        if compounds:
            compound = compounds[0] # Take the first result
            inchi = compound.inchi
            return {"name": name, "cid": compound.cid, "inchi": inchi}
        # if not found, try alternative names
        else:
            return try_alternative_name(name)

    
    parse_chemical = sc(parse_chemical)
    results = []
    for name in tqdm(chemical_names, desc="Parsing chemicals", unit="chemical"):
        res = parse_chemical(name)
        if res is not None:
            results.append(res)
    
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    logging.info(f"Saved results to {output_path}")
    
    return df


