import pandas as pd
import logging
import pathlib
import time 
import toxindex.utils.chemprop as chemprop
import itertools
import toxindex.utils.simplecache as simplecache
from tqdm import tqdm

def predict_chemicals(input_path, output_path):
    """
    Predict chemical properties using chemprop and save results.
    
    Args:
        input_path (str or pathlib.Path): Path to the input CSV file containing chemicals with 'inchi' column
        output_path (str or pathlib.Path): Path to save the output parquet file with predictions
    """
    # Convert paths to pathlib.Path objects if they're not already
    input_path = pathlib.Path(input_path)
    output_path = pathlib.Path(output_path)
    
    # Read input data
    indf = pd.read_csv(input_path)
    
    # Setup cache for predictions
    cachedir = pathlib.Path('cache') / 'function_cache' / 'chemprop_predictions_cache'
    cachedir.mkdir(parents=True, exist_ok=True)
    pred = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)

    def safe_predict(inchi):
        try:
            return pred(inchi) if inchi else None
        except Exception as e:
            logging.error(f"Error predicting for inchi {inchi}: {e}")
            return []
    
    # Make predictions
    predictions = [safe_predict(inchi) for inchi in tqdm(indf['inchi'])]
    predictions = list(itertools.chain.from_iterable(predictions))
    pdf = pd.DataFrame(predictions)
    
    # Extract property info
    pdf['property_title'] = pdf['property'].apply(lambda x: str(x.get('title', ''))).astype(str)
    pdf['property_source'] = pdf['property'].apply(lambda x: str(x.get('source', ''))).astype(str)
    pdf['property_categories'] = pdf['property'].apply(lambda x: str(x.get('categories', ''))).astype(str)
    pdf['property_metadata'] = pdf['property'].apply(lambda x: str(x.get('metadata', ''))).astype(str)
    pdf['property'] = pdf['property'].astype(str)
    
    # Merge with categorization data
    resdf = pd.merge(pdf, indf, on='inchi', how='left')
    resdf.to_parquet(output_path)
    
    logger.info(f"Saved results to {output_path}")
    
    return resdf
