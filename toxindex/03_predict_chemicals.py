import pandas as pd
import logging
import pathlib
import time 
import toxindex.utils.chemprop as chemprop
import itertools
import toxindex.utils.simplecache as simplecache
from tqdm import tqdm

# --- DEPENDENCIES ---
indf = pd.read_csv(pathlib.Path('cache/categorize_chemicals/classified_chemicals.csv'))
catdf = pd.read_csv(pathlib.Path('cache/categorize_chemicals/classified_chemicals.csv'))

# --- Configuration ---
outdir = pathlib.Path('cache/predict_chemicals') # Use 'outdir', new path
indf = pd.read_csv(pathlib.Path('cache/categorize_chemicals/classified_chemicals.csv'))

# --- Setup Output Directory and Logging ---
outdir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=outdir / 'log.txt', # New log path
    filemode='w',                # Overwrite log
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


cachedir = outdir / 'chemprop_predictions_cache'
cachedir.mkdir(parents=True, exist_ok=True)
pred = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)
predictions = [pred(inchi) for inchi in tqdm(indf['inchi'])]
predictions = list(itertools.chain.from_iterable(predictions))
pdf = pd.DataFrame(predictions)

# Extract property info
pdf['property_title'] = pdf['property'].apply(lambda x: str(x.get('title', ''))).astype(str)
pdf['property_source'] = pdf['property'].apply(lambda x: str(x.get('source', ''))).astype(str)
pdf['property_categories'] = pdf['property'].apply(lambda x: str(x.get('categories', ''))).astype(str)
pdf['property_metadata'] = pdf['property'].apply(lambda x: str(x.get('metadata', ''))).astype(str)
pdf['property'] = pdf['property'].astype(str)

# get the categorization 
resdf = pd.merge(pdf, catdf, on='inchi', how='left')
resdf.to_parquet(outdir / 'chemprop_predictions.parquet')

logger.info(f"Saved results to {outdir}")


