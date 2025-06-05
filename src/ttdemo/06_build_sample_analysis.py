import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image, ImageOps # ImageOps might be needed for inversion later if direct numpy fails
import io
import warnings # To ignore tight_layout warning

# Attempt to import Cairo drawer
try:
    from rdkit.Chem.Draw import MolDraw2DCairo
except ImportError:
    print("ERROR: MolDraw2DCairo not available. Please install 'cairocffi'.")
    exit()

# --- Configuration ---
CACHE_DIR = pathlib.Path('cache')
RESOURCES_DIR = CACHE_DIR / 'resources'
PREDICT_DIR = CACHE_DIR / 'predict_chemicals'
OUTPUT_DIR = CACHE_DIR / 'build_sample_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Load Data ---
gemini_props = pd.read_csv(RESOURCES_DIR / 'gemini-interesting-properties.txt', header=None)[0].tolist()
gemini_chems = pd.read_csv(RESOURCES_DIR / 'gemini-interesting-chemicals', header=None)[0].tolist()
current_date_str = pd.Timestamp.now(tz='America/New_York').strftime('%Y-%m-%d %H:%M:%S %Z')
print(f"Analysis run at: {current_date_str}")

df = pd.read_parquet(PREDICT_DIR / 'chemprop_predictions.parquet')

# --- Filter Data ---
df = df[df['property_title'].isin(gemini_props) & df['name'].isin(gemini_chems)].copy()


# --- Data Preparation ---
# Assign property numbers based on the original gemini_props order
prop_to_num = {prop: i + 1 for i, prop in enumerate(gemini_props)}
num_to_prop = {v: k for k, v in prop_to_num.items()} # For potential later use
df['property_num'] = df['property_title'].map(prop_to_num)

# Calculate mean activity per chemical using ORIGINAL values
chem_mean_activity = df.groupby('name')['value'].mean().sort_values().reset_index()
chem_order = chem_mean_activity['name'].tolist()
n_chems = len(chem_order)

# Pivot table with ORIGINAL values for the main heatmap
pivot_df = df.pivot_table(index='name', columns='property_title', values='value', aggfunc='first')
pivot_df = pivot_df.reindex(chem_order) # Reorder rows

# --- START: Ensure Correct Column Order (p1 to p6) ---

# 1. Create mapping from original property titles to p1, p2, etc.
prop_mapping = {old_col: f'p{prop_to_num[old_col]}' for old_col in pivot_df.columns}

# 2. Rename the columns
pivot_df = pivot_df.rename(columns=prop_mapping)

# 3. Define the desired column order
num_properties = len(gemini_props) # Should be 6 in this case
desired_column_order = [f'p{i+1}' for i in range(num_properties)]

# 4. Reindex the DataFrame columns using the desired order
# Check if all desired columns exist before reindexing to avoid errors
existing_cols_in_order = [col for col in desired_column_order if col in pivot_df.columns]
pivot_df = pivot_df[existing_cols_in_order]


# Prepare data for the mean activity heatmap column
mean_heatmap_df = chem_mean_activity.set_index('name')[['value']].reindex(chem_order)
mean_heatmap_df = mean_heatmap_df.rename(columns={'value': 'Mean Activity'})

# --- Plotting Setup ---
plt.rcParams.update({
    'figure.facecolor': 'none', 'axes.facecolor': 'none', 'savefig.facecolor': 'none',
    'text.color': 'white', 'axes.labelcolor': 'white', 'axes.edgecolor': 'white',
    'xtick.color': 'white', 'ytick.color': 'white',
})
sns.set(style="dark")

# --- Create Figure and Axes ---
fig = plt.figure(figsize=(22, 10))
# Wider main heatmap relative width in gridspec
gs = fig.add_gridspec(1, 3, width_ratios=[len(pivot_df.columns) + 2, 1, 0.5], wspace=0.1)

ax_main_heat = fig.add_subplot(gs[0, 0])
# REMOVED sharey=ax_main_heat
ax_mean_heat = fig.add_subplot(gs[0, 1])
cbar_ax = fig.add_subplot(gs[0, 2])

# --- Annotation Color Choice ---
annotation_color = 'black'
annotation_kws = {"size": 9, "color": annotation_color, "weight": "bold"}

# --- Main Heatmap ---
heatmap_main = sns.heatmap(
    pivot_df, cmap="viridis", linewidths=0.8, linecolor='black',
    annot=True, fmt=".2f", cbar=True, cbar_ax=cbar_ax,
    cbar_kws={'label': 'Activity Value'}, ax=ax_main_heat, annot_kws=annotation_kws
)
ax_main_heat.set_ylabel("")
ax_main_heat.set_xlabel("Properties", fontsize=12, color='white')
ax_main_heat.set_title("Chemical Activity Values", fontsize=16, color='white')

# Ensure Y tick labels (chemical names) are visible and WHITE
ax_main_heat.tick_params(axis='y', labelcolor='white', labelsize=10) # Reduced size slightly
# Force labels white after plotting
for label in ax_main_heat.get_yticklabels():
    label.set_color('white')

# Relabel X-axis
ax_main_heat.set_xticks(np.arange(len(pivot_df.columns)) + 0.5)
ax_main_heat.set_xticklabels(pivot_df.columns, fontsize=10, color='white', rotation=0)

# Style color bar
cbar = ax_main_heat.collections[0].colorbar
cbar.ax.yaxis.set_tick_params(color='white', labelsize=10)
cbar.outline.set_edgecolor('white')
cbar.set_label('Activity Value', color='white', fontsize=12)

# --- Mean Activity Heatmap ---
heatmap_mean = sns.heatmap(
    mean_heatmap_df, cmap="viridis", linewidths=0.8, linecolor='black',
    annot=True, fmt=".2f", cbar=False, ax=ax_mean_heat, annot_kws=annotation_kws
)
ax_mean_heat.set_yticklabels([]) # Explicitly remove Y labels here
ax_mean_heat.set_ylabel("")
ax_mean_heat.set_xlabel("")
ax_mean_heat.set_title("Mean\nActivity", fontsize=12, color='white')

# Set X tick labels for mean heatmap
ax_mean_heat.set_xticks([0.5])
ax_mean_heat.set_xticklabels([mean_heatmap_df.columns[0]], rotation=0, color='white', fontsize=10)
ax_mean_heat.tick_params(bottom=False, labelbottom=True)

# --- Manual Y-axis Alignment (since sharey was removed) ---
# Do this AFTER drawing both heatmaps
fig.canvas.draw() # Ensure layout is calculated
ax_mean_heat.set_ylim(ax_main_heat.get_ylim())


# --- Save Plot ---
# Rely mainly on tight_layout, adjust rect slightly if needed
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    # Give more space on left with rect=[0.1, ...]
    fig.tight_layout(rect=[0.1, 0, 0.95, 1])
# REMOVED subplots_adjust

plt.savefig(OUTPUT_DIR / 'sample_chemical_heatmap_unnormalized.png', dpi=300, facecolor='none', edgecolor='none')
plt.close(fig)
print(f"Unnormalized heatmap saved to {OUTPUT_DIR / 'sample_chemical_heatmap_unnormalized.png'}")


# --- Draw Chemical Structures (Cairo + Post-Inversion) ---
def draw_chemicals_cairo_then_invert(name_inchi_pairs, outfile, mols_per_row=3, img_size=(300, 300)):
    mols = []
    valid_names = []
    atomic_nums = set()
    for name, inchi in name_inchi_pairs:
        mol = Chem.MolFromInchi(inchi)
        # NO error handling
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        mol = Chem.RemoveHs(mol)
        mols.append(mol)
        valid_names.append(name)
        for atom in mol.GetAtoms():
            atomic_nums.add(atom.GetAtomicNum())

    n_mols = len(mols)
    if n_mols == 0:
        print("No molecules to draw.")
        return

    n_cols = min(mols_per_row, n_mols)
    n_rows = (n_mols + n_cols - 1) // n_cols
    panel_width, panel_height = img_size

    # Initialize Cairo drawer
    drawer = MolDraw2DCairo(n_cols * panel_width, n_rows * panel_height, panel_width, panel_height)
    opts = drawer.drawOptions()

    # Configure Cairo to draw WHITE on TRANSPARENT (as before)
    opts.setBackgroundColour((0.0, 0.0, 0.0, 0.0))
    white = (1.0, 1.0, 1.0)
    opts.useDefaultAtomPalette = False
    white_palette = {num: white for num in atomic_nums}
    white_palette[0] = white
    opts.atomPalette = white_palette
    opts.legendColour = white

    # Draw molecules
    drawer.DrawMolecules(mols, legends=valid_names)
    drawer.FinishDrawing()

    # Get image data from Cairo
    png_data = drawer.GetDrawingText()
    img_from_cairo = Image.open(io.BytesIO(png_data))

    # --- POST-PROCESSING: Ensure foreground is white ---
    img_rgba = img_from_cairo.convert('RGBA')
    data = np.array(img_rgba)

    # Mask for non-transparent pixels (alpha > threshold)
    alpha_channel = data[:, :, 3]
    opaque_mask = alpha_channel > 10 # Pixels belonging to structure/text

    # Set RGB of opaque pixels to WHITE (255, 255, 255)
    data[opaque_mask, 0:3] = 255

    # Create final image from modified data
    img_final_white = Image.fromarray(data)
    # --------------------------------------------------

    # Save the final (explicitly whitened) image
    img_final_white.save(outfile, format="PNG")
    print(f"Chemical structures (white on transparent) saved to {outfile}")


chemdf = df[['name', 'inchi']].drop_duplicates()
name_inchi_pairs = list(zip(chemdf['name'], chemdf['inchi']))
# Call the drawing function
draw_chemicals_cairo_then_invert(name_inchi_pairs, OUTPUT_DIR / 'chemical_structures_white.png')


# --- Create and Save Summary Table ---
summary_df = chem_mean_activity.rename(columns={'name': 'Chemical', 'value': 'Mean Activity'})
summary_df['Mean Activity'] = summary_df['Mean Activity'].round(4)
summary_df = summary_df[['Chemical', 'Mean Activity']]
summary_df.to_csv(OUTPUT_DIR / 'chemical_activity_summary.csv', index=False)
print(f"Summary saved to {OUTPUT_DIR / 'chemical_activity_summary.csv'}")

print("\nAnalysis complete!")