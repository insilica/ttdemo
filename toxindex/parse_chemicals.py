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
    def parse_chemical(name):
        logging.info(f"Parsing chemical: {name}")
        time.sleep(0.33)
        compounds = pcp.get_compounds(name, 'name')
        if not compounds:
            return None
        compound = compounds[0] # Take the first result
        inchi = compound.inchi
        return {"name": name, "cid": compound.cid, "inchi": inchi}
    
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
