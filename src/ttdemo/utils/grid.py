from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import matplotlib.pyplot as plt

def create_grid_image(molecules_df, value_col='value', n_per_row=4, output_path=None, title=None):
    """Create a grid image of molecules from InChI strings."""
    mols = []
    legends = []
    
    for _, row in molecules_df.iterrows():
        name = row['name']
        inchi = row['inchi']
        value = row[value_col]
        
        # Create molecule from InChI
        mol = Chem.MolFromInchi(inchi)
        
        # Remove hydrogens and compute 2D coordinates
        mol = Chem.RemoveAllHs(mol)
        AllChem.Compute2DCoords(mol)
        mol = Draw.PrepareMolForDrawing(mol)
        mols.append(mol)
        
        # Create legend
        short_name = name if len(name) < 20 else name[:17] + "..."
        legend = f"{short_name}\nValue: {value:.2f}"
        legends.append(legend)
    
    # Create the grid image with basic settings
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=n_per_row,
        subImgSize=(250, 250),
        legends=legends,
        useSVG=False
    )
    
    # Add title and save
    if title:
        plt.figure(figsize=(10, 10))
        plt.imshow(img)
        plt.title(title, fontsize=16)
        plt.axis('off')
        if output_path:
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        return img
    
    if output_path and not title:
        img.save(output_path)
    
    return img