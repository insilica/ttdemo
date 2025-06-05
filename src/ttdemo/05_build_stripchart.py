import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

def build_stripchart(ringdf, outfile):
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
    # Set the style to dark background
    plt.style.use('dark_background')
    
    # Define color mapping for categories (matching the heatmap)
    category_colors = {
        '1 Ring Aromatic': '#FFFED0', '2 Ring Aromatic': '#FBDA80', 
        '3 Ring Aromatic': '#F7B659', '4 Ring Aromatic': '#EE6033',
        '5 Ring Aromatic': '#D53D23', '6+ Ring Aromatic': '#781A26',
    }
    
    # Create a more predictable category order
    category_order = [
        '1 Ring Aromatic', '2 Ring Aromatic',
        '3 Ring Aromatic', '4 Ring Aromatic', '5 Ring Aromatic', '6+ Ring Aromatic'
    ]
    
    # Create the figure
    plt.figure(figsize=(12, 7))
    
    # Create the stripplot
    ax = sns.stripplot(x='classification', y='value', data=ringdf, 
                      palette=category_colors, size=8, jitter=True, alpha=0.8,
                      order=category_order)
    
    # Add horizontal lines for means
    means = ringdf.groupby('classification')['value'].mean()
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
    plt.savefig(outfile, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Strip chart saved to {outfile}")

gemini_props = pd.read_csv(pathlib.Path('cache/resources/gemini-properties.txt'), header=None)
df = pd.read_parquet(pathlib.Path('cache/predict_chemicals/chemprop_predictions.parquet'))
df = df.merge(gemini_props, left_on='property_title', right_on=0, how='inner')
stripdf = ringdf.groupby(['classification','name'])['value'].mean().reset_index()
build_stripchart(stripdf, outdir / 'ring_stripchart.png')