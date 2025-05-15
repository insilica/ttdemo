import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

import logging

logger = logging.getLogger(__name__)

def build_stripchart(input_path, output_path):
    """
    Create a strip chart showing mean property values by chemical classification.
    
    Parameters:
    -----------
    ringdf : pandas DataFrame
        DataFrame containing the chemical data with columns:
        - 'classification': Chemical class (e.g. '1 Ring System', 'C9 paraffin')
        - 'value': Numerical property values
    outfile : str or pathlib.Path
        Path to save the output plot
        
    Returns:
    --------
    None
    """
    # Convert paths to pathlib.Path objects if they're not already
    input_path = pathlib.Path(input_path)
    output_path = pathlib.Path(output_path)

    # Set the style to dark background
    plt.style.use('dark_background')

    # gemini_props = pd.read_csv(pathlib.Path('cache/resources/gemini-properties.txt'), header=None)
    df_allfeat = pd.read_parquet(input_path)

    feature_path = input_path.parent / 'matched_properties.txt'
    feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    feature_names = sorted(list(feature_names))
    df_allfeat['is_in_lookup'] = df_allfeat['property_title'].isin(feature_names)
    df = df_allfeat[df_allfeat['is_in_lookup']]

    classdf = pd.read_csv(input_path.parent / 'classified_chemicals.csv')
    if 'classification' not in df.columns:
        df = df.merge(classdf, on='inchi', how='left')

    if 'name' not in df.columns:
        # print(df.columns)
        df['name'] = df['name_x']
    # df = df.merge(gemini_props, left_on='property_title', right_on=0, how='inner')
    stripdf = df.groupby(['classification','name'])['value'].mean().reset_index()
    
    # Define color mapping for categories (matching the heatmap)
    # category_colors = {
    #     '1 Ring Aromatic': '#FFFED0', '2 Ring Aromatic': '#FBDA80', 
    #     '3 Ring Aromatic': '#F7B659', '4 Ring Aromatic': '#EE6033',
    #     '5 Ring Aromatic': '#D53D23', '6+ Ring Aromatic': '#781A26',
    # }
    category_colors = {'Nephrotoxic': '#FFFED0', 'Non-Nephrotoxic': '#FBDA80'}

    # # Create a more predictable category order
    # category_order = [
    #     '1 Ring Aromatic', '2 Ring Aromatic',
    #     '3 Ring Aromatic', '4 Ring Aromatic', '5 Ring Aromatic', '6+ Ring Aromatic'
    # ]
    category_order = ['Nephrotoxic','Non-Nephrotoxic']
    
    # Create the figure
    plt.figure(figsize=(12, 7))
    
    # Create the stripplot
    ax = sns.stripplot(x='classification', y='value', data=stripdf, 
                      palette=category_colors, size=8, jitter=True, alpha=0.8,
                      order=category_order)
    
    # Add horizontal lines for means
    means = stripdf.groupby('classification')['value'].mean()
    for cat in category_order:
        if cat in means:
            plt.hlines(y=means[cat], 
                      xmin=ax.get_xticks()[category_order.index(cat)] - 0.4,
                      xmax=ax.get_xticks()[category_order.index(cat)] + 0.4, 
                      colors='magenta', linewidth=2)
    
    # Customize the plot
    plt.ylim(0, 0.7)  # Set y-axis limits
    plt.ylabel('Mean property value', fontsize=12)
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Add a descriptive title
    plt.title('Top Properties by Activity Level', fontsize=16, pad=20)
    
    # Save and close the plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Strip chart saved to {output_path}")

if __name__ == "__main__":
    gemini_props = pd.read_csv(pathlib.Path('cache/resources/gemini-properties.txt'), header=None)
    df = pd.read_parquet(pathlib.Path('cache/predict_chemicals/chemprop_predictions.parquet'))
    df = df.merge(gemini_props, left_on='property_title', right_on=0, how='inner') #filter 
    stripdf = ringdf.groupby(['classification','name'])['value'].mean().reset_index()
    build_stripchart(stripdf, outdir / 'ring_stripchart.png')