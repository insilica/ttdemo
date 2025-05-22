import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import yaml
import logging
from scipy.stats import ttest_ind

logger = logging.getLogger(__name__)

def build_stripchart(input_path, output_path, agg_func, feature_selection_method):
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

    # feature_path = input_path.parent / 'matched_properties.txt'
    feature_path = input_path.parent / 'selected_properties' / f'{feature_selection_method}_selected_properties.txt'
    feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    feature_names = sorted(list(feature_names))
    df_allfeat['is_in_lookup'] = df_allfeat['property_title'].isin(feature_names)
    df = df_allfeat[df_allfeat['is_in_lookup']]

    classdf = pd.read_csv(input_path.parent / 'classified_chemicals.csv')
    if 'classification' not in df.columns or df['classification'].isna().any():
        df = df.drop(columns=['classification'], errors='ignore')  # drop to avoid _x/_y
        df = df[df['inchi'].isin(classdf['inchi'])] #filter inchi with classified label
        df = df.merge(classdf, on='inchi', how='left')
        if 'classification' not in df.columns:
            raise ValueError("Merge failed: 'classification' column is missing after merging with 'classified_chemicals.csv'.") 

    if 'name' not in df.columns:
        # print(df.columns)
        df['name'] = df['name_x']
    # df = df.merge(gemini_props, left_on='property_title', right_on=0, how='inner')
    stripdf = df.groupby(['classification', 'name'])['value'].agg(agg_func).reset_index()
    

    # Get unique classifications
    unique_classes = sorted(df['classification'].unique())  # sorted for consistency

    # Load colors from YAML
    with open("cache/resources/colormap.yaml", "r") as f:
        config = yaml.safe_load(f)

    colors_hex = config['colors']

    # Safely assign only as many colors as needed
    category_colors = dict(zip(unique_classes, colors_hex[:len(unique_classes)]))
    category_order = unique_classes

    # Define color mapping for categories (matching the heatmap)
    # category_colors = {
    #     '1 Ring Aromatic': '#FFFED0', '2 Ring Aromatic': '#FBDA80', 
    #     '3 Ring Aromatic': '#F7B659', '4 Ring Aromatic': '#EE6033',
    #     '5 Ring Aromatic': '#D53D23', '6+ Ring Aromatic': '#781A26',
    # }
    # category_colors = {'DNT': '#FFFED0', 'Non-DNT': '#FBDA80'}
    
    # category_colors = {'hepatotoxic': '#FFFED0', 'nontoxic': '#FBDA80'}

    # # Create a more predictable category order
    # category_order = [
    #     '1 Ring Aromatic', '2 Ring Aromatic',
    #     '3 Ring Aromatic', '4 Ring Aromatic', '5 Ring Aromatic', '6+ Ring Aromatic'
    # ]
    # category_order = ['DNT','Non-DNT']
    # category_order = ['Nephrotoxic','Non-nephrotoxic']
    # category_order = ['hepatotoxic','nontoxic']
    
    # Create the figure
    plt.figure(figsize=(12, 7))
    
    # Create the stripplot
    # ax = sns.stripplot(x='classification', y='value', data=stripdf, 
    #                   palette=category_colors, size=8, jitter=True, alpha=0.8,
    #                   order=category_order)
    
    ax = sns.stripplot(
    x='classification',
    y='value',
    hue='classification',  # explicitly define hue
    data=stripdf,
    palette=category_colors,
    size=8,
    jitter=True,
    alpha=0.8,
    order=category_order,
    legend=False  # avoid duplicate legend
    )

    classes = df['classification'].unique()
    if len(classes) != 2:
        raise ValueError("T-test requires exactly two groups for comparison.")

    # Split the groups
    group_a = stripdf[stripdf['classification'] == classes[0]]['value']
    group_b = stripdf[stripdf['classification'] == classes[1]]['value']

    # Perform independent t-test (assume unequal variance just to be safe)
    t_stat, p_value = ttest_ind(group_a, group_b, equal_var=False)

    print(f"T-statistic: {t_stat:.4f}")
    print(f"P-value: {p_value:.4e}")

    # Check for significance (e.g., alpha = 0.05)
    if p_value < 0.05:
        print("Result: Significant difference between groups")
    else:
        print("Result: No significant difference between groups")
    
    # Add horizontal lines for means
    if agg_func == 'mean':
        stats = stripdf.groupby('classification')['value'].mean()
    elif agg_func == 'median':
        stats = stripdf.groupby('classification')['value'].median()
    else:
        raise ValueError("agg_func must be 'mean' or 'median'")

    for cat in category_order:
        if cat in stats:
            plt.hlines(y=stats[cat], 
                    xmin=ax.get_xticks()[category_order.index(cat)] - 0.4,
                    xmax=ax.get_xticks()[category_order.index(cat)] + 0.4, 
                    colors='magenta', linewidth=2)
    
    # Customize the plot
    # Get current y-limits from the data
    ymin, ymax = plt.ylim()
    plt.ylim(0, ymax * 1.1)  # Add 10% space above the max

    plt.ylabel(f"{agg_func.capitalize()} property value", fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=20)
    plt.tight_layout()
    
    # Add a descriptive title
    plt.title('Top Properties by Activity Level', fontsize=24, pad=20)
    
    # Save and close the plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    
    print(f"Strip chart saved to {output_path}")

if __name__ == "__main__":
    gemini_props = pd.read_csv(pathlib.Path('cache/resources/gemini-properties.txt'), header=None)
    df = pd.read_parquet(pathlib.Path('cache/predict_chemicals/chemprop_predictions.parquet'))
    df = df.merge(gemini_props, left_on='property_title', right_on=0, how='inner') #filter 
    stripdf = ringdf.groupby(['classification','name'])['value'].mean().reset_index()
    build_stripchart(stripdf, outdir / 'ring_stripchart.png', stats)