import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
from sklearn.preprocessing import MinMaxScaler
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import logging

logger = logging.getLogger(__name__)

def build_heatmap(input_path, output_path):
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

    # filter feature
    feature_path = input_path.parent / 'claude_relevant_properties.txt'
    feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    feature_names = sorted(list(feature_names))
    df_allfeat['is_in_lookup'] = df_allfeat['property_title'].isin(feature_names)
    df = df_allfeat[df_allfeat['is_in_lookup']]
    
    classdf = pd.read_csv(input_path.parent / 'classified_chemicals.csv')
    if 'classification' not in df.columns:
        df = df.merge(classdf, on='inchi', how='left')

    if 'name' not in df.columns:
        print(df.columns)
        df['name'] = df['name_x']
        
    # Input for heatmap
    pdf = df[['name', 'property_token', 'value', 'classification']]

    # Filter data for ring compounds
    # pdf = pdf[pdf['classification'].str.contains("Ring")]
    
    # Generate the heatmap
    _generate_heatmap(pdf, output_path)
    
    # Create a csv of the top 10 most activated properties
    top_props_path = output_path.parent / 'top_props.csv'
    top10df = df[['name', 'property_title', 'property_source', 'property_metadata', 'value', 'classification']]
    top10df = top10df[~top10df['classification'].str.contains("Paraffin")]
    
    top10_props = top10df\
        .groupby(['property_title', 'property_source', 'property_metadata'])['value'].mean()\
        .reset_index().sort_values('value', ascending=False)\
        .head(10)
    
    top10_props.to_csv(top_props_path, index=False)
    logger.info(f"Saved top properties to {top_props_path}")
    
    return output_path


def _generate_heatmap(pdf, outfile):
    """
    Internal function to generate the actual heatmap visualization.
    
    Args:
        pdf (pandas.DataFrame): DataFrame containing the chemical prediction data
        outfile (pathlib.Path): Path to save the output heatmap image
    """
    # Pivot and normalize
    pivot_df = pdf.pivot_table(index='name', columns='property_token', values='value', aggfunc='first')
    norm_values = pd.DataFrame(
        MinMaxScaler().fit_transform(pivot_df.fillna(0)), 
        index=pivot_df.index, 
        columns=pivot_df.columns
    )

    # Direct color mapping (no need for intermediate classification)
    # category_colors = {
    #     '1 Ring Aromatic': '#FFFED0', '2 Ring Aromatic': '#FBDA80', 
    #     '3 Ring Aromatic': '#F7B659', '4 Ring Aromatic': '#EE6033',
    #     '5 Ring Aromatic': '#D53D23', '6+ Ring Aromatic': '#781A26',
    # }
    category_colors = {'hepatotoxic': '#FFFED0', 'nontoxic': '#FBDA80'}

    # Apply classification and get row colors
    class_series = pdf.drop_duplicates('name').set_index('name')['classification']
    norm_values['substance_class'] = class_series
    row_colors = norm_values['substance_class'].map(category_colors)
    norm_values = norm_values.drop(columns='substance_class')

    # Create a wider figure to accommodate the colorbar
    plt.figure(figsize=(18, 9))  # Increased width from 16 to 18
    
    # Plot heatmap with adjusted layout - keep column clustering
    sns.set(style="white")

    # if the map shows black, reduce linewidths.
    g = sns.clustermap(
        norm_values, cmap="viridis", row_colors=row_colors,
        xticklabels=False, yticklabels=False, 
        # linewidths=0.001,
        linecolor='black', col_cluster=True, row_cluster=True,
        figsize=(18, 9),  # Wider figure
        cbar_pos=(0.91, 0.3, 0.02, 0.4),  # More to the left and higher up
        dendrogram_ratio=(0.1, 0.05),  # Make column dendrogram shorter (was 0.2 by default)
        tree_kws={'linewidths': 0.5}  # Thinner dendrogram lines
    )
    
    # Hide column dendrogram but still keep clustering
    g.ax_col_dendrogram.set_visible(False)
    
    # Adjust the main heatmap's position to be more centered and shifted right
    heatmap_pos = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0,  # No left shift, keep original x position
                               heatmap_pos.width * 0.88, heatmap_pos.height])  # Slightly narrower
    
    # Adjust the row dendrogram position as well
    row_pos = g.ax_row_dendrogram.get_position()
    g.ax_row_dendrogram.set_position([row_pos.x0, row_pos.y0,  # Keep original x position
                                     row_pos.width, row_pos.height])
    
    # Adjust row colors panel position
    if hasattr(g, 'ax_row_colors'):
        row_colors_pos = g.ax_row_colors.get_position()
        g.ax_row_colors.set_position([row_colors_pos.x0, row_colors_pos.y0,  # Keep original x position
                                     row_colors_pos.width, row_colors_pos.height])

    # Create a separate axis for the class legend on the left (moved further left)
    left_legend_ax = plt.axes([0.005, 0.7, 0.05, 0.2])  # [x, y, width, height] - moved x from 0.01 to 0.005
    left_legend_ax.axis('off')  # Hide the axis

    # Add class legend to the left axis
    class_legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in category_colors.values()]
    left_legend_ax.legend(
        class_legend_handles, 
        category_colors.keys(), 
        loc='center', 
        frameon=False, 
        title='Class', 
        title_fontsize=10,
        labelcolor='white',
        prop={'size': 8}
    )

    # Adjust the colorbar (scale legend) that's now on the far right
    cbar = g.ax_cbar
    cbar.set_ylabel('Scale', color='white', fontsize=10)
    cbar.tick_params(colors='white', labelsize=8)
    
    # Add axis labels but no title
    g.ax_heatmap.set_xlabel("Estimated property values", color='white')
    g.ax_heatmap.set_ylabel("Substance", color='white')
    
    plt.savefig(outfile, dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved heatmap to {outfile}")