import pandas as pd
import logging
import pathlib
import time
import toxindex.utils.chemprop as chemprop
import itertools
import toxindex.utils.simplecache as simplecache
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def predict_chemicals(input_csv, outdir, max_workers=30):
    # --- Load Data ---
    indf = pd.read_csv(input_csv)
    catdf = pd.read_csv(input_csv)

    # --- Setup Output Directory and Logging ---
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=outdir / "log.txt",
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)

    # --- Setup Cache ---
    cachedir = outdir / "chemprop_predictions_cache"
    cachedir.mkdir(parents=True, exist_ok=True)
    pred = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)

    # --- Parallel Prediction ---
    predictions = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(pred, row["inchi"]): row["inchi"]
            for _, row in indf.iterrows()
        }
        for future in tqdm(as_completed(futures), total=len(futures), desc="Predicting"):
            inchi = futures[future]
            try:
                result = future.result()
                predictions.extend(result)
            except Exception as e:
                logger.error(f"Failed to predict {inchi}: {e}")

    pdf = pd.DataFrame(predictions)

    # --- Extract property info ---
    pdf["property_title"] = pdf["property"].apply(lambda x: str(x.get("title", "")))
    pdf["property_source"] = pdf["property"].apply(lambda x: str(x.get("source", "")))
    pdf["property_categories"] = pdf["property"].apply(lambda x: str(x.get("categories", "")))
    pdf["property_metadata"] = pdf["property"].apply(lambda x: str(x.get("metadata", "")))
    pdf["property"] = pdf["property"].astype(str)

    # --- Merge & Save ---
    logger.info(f"Saving results to {outdir}")
    resdf = pd.merge(pdf, catdf, on="inchi", how="left")
    resdf.to_parquet(outdir / "chemprop_predictions.parquet")

    logger.info(f"Saved results to {outdir}")
