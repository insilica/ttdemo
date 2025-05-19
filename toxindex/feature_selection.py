import logging
import pathlib
import pandas as pd

logger = logging.getLogger(__name__)

def select_feature(input_path, output_path):
    """
    Build a heatmap visualization from chemical prediction data.
    
    Args:
        input_path (str or pathlib.Path): Path to the parquet file containing prediction data
        output_path (str or pathlib.Path): Path to save the output heatmap image and related files
    """
    # Convert paths to pathlib.Path objects if they're not already
    input_path = pathlib.Path(input_path)
    output_path = pathlib.Path(output_path)
    
    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
# tell pandas to never wrap, just keep going wider
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 100)
    
    logger.info(f"Building heatmap from {input_path}")
    
    # Load prediction data
    df_allfeat = pd.read_parquet(input_path)

    classdf = pd.read_csv(input_path.parent / 'classified_chemicals.csv')
    # if 'classification' not in df.columns:
    #     df = df.merge(classdf, on='inchi', how='left')

    if 'classification' not in df.columns or df['classification'].isna().any():
    # Merge classification info
        print('reading classification')
        df = df.drop(columns=['classification'], errors='ignore')  # drop to avoid _x/_y
        df = df.merge(classdf[['inchi', 'classification']], on='inchi', how='left')

    if 'name' not in df.columns:
        # print(df.columns)
        df['name'] = df['name_x']
        
    # Input for heatmap
    pdf = df[['name', 'property_token', 'value', 'classification']]

    pivot_df = pdf.pivot_table(index='name', columns='property_token', values='value', aggfunc='first')
    norm_values = pd.DataFrame(
        MinMaxScaler().fit_transform(pivot_df.fillna(0)), 
        index=pivot_df.index, 
        columns=pivot_df.columns
    )

    class_series = pdf.drop_duplicates('name').set_index('name')['classification']
    norm_values['substance_class'] = class_series








    #%% feature selection



    ###


    # filter feature
    feature_path = input_path.parent / 'matched_properties.txt'
    feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    feature_names = sorted(list(feature_names))
    df_allfeat['is_in_lookup'] = df_allfeat['property_title'].isin(feature_names)
    df = df_allfeat[df_allfeat['is_in_lookup']]

    


if __name__ == "__main__":
    project = 'hepatotoxic'

    cachedir = pathlib.Path('cache')
    cachedir.mkdir(exist_ok=True)

    outdir = cachedir / 'projects' / project / 'heatmap_dir'
    outdir.mkdir(exist_ok=True)

    input_path=cachedir / 'projects' / project / 'predictions.parquet'
    output_path=outdir / 'heatmap.png'
    select_feature(input_path, output_path)