import rapidfuzz
import numpy as np
import pathlib
import json
import ttdemo.build_stripchart as build_stripchart
import ttdemo.build_heatmap as build_heatmap
import ttdemo.parse_chemicals as parse_chemicals
import ttdemo.predict_chemicals as predict_chemicals
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

def contains_sulfonate(inchi: str) -> bool:
    """
    Detects if a molecule contains a sulfonate group using SMARTS.
    Works for -SO3H, -SO3-, and sulfonate esters.
    """
    mol = Chem.MolFromInchi(inchi)
    sulfonate_smarts = Chem.MolFromSmarts("S(=O)(=O)[O-,$([OH]),$([O][CX4])]")
    return mol.HasSubstructMatch(sulfonate_smarts)

def get_perfluorinated_chain_length(inchi: str) -> int:
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return 0

    # Identify all carbon atoms that are only bonded to C or F
    eligible_carbons = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # carbon
            continue
        neighbors = atom.GetNeighbors()
        if all(nbr.GetAtomicNum() in (6, 9) for nbr in neighbors):  # C or F only
            eligible_carbons.add(atom.GetIdx())

    # Build the longest continuous chain from eligible carbon atoms
    from collections import deque

    def bfs(start_idx):
        visited = set()
        queue = deque([(start_idx, 1)])
        max_len = 1
        while queue:
            idx, length = queue.popleft()
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in eligible_carbons and nbr_idx not in visited:
                    queue.append((nbr_idx, length + 1))
                    max_len = max(max_len, length + 1)
        return max_len

    max_chain = 0
    for start in eligible_carbons:
        chain_len = bfs(start)
        max_chain = max(max_chain, chain_len)

    return max_chain

cachedir = pathlib.Path('cache') / 'notebooks' / 'pfas_analysis'
cachedir.mkdir(parents=True,exist_ok=True)

# get chemicals
ipath = cachedir / "chemicals.txt"
opath = cachedir / "parsed_chemicals.csv"
parse_chemicals.parse_chemicals(input_path=ipath,output_path=opath)

# find chain lengths
chemicals = pd.read_csv(cachedir / "parsed_chemicals.csv")
chemicals['chain_length'] = chemicals['inchi'].apply(get_perfluorinated_chain_length)
chemicals['sulfonate'] = chemicals['inchi'].apply(contains_sulfonate)

# predict properties
ipath = cachedir / "parsed_chemicals.csv"
opath = cachedir / "predictions.parquet"
predict_chemicals.predict_chemicals(input_path=ipath,output_path=opath)

# get pfas and hepatotoxicity  predictions
pfas_predictions = pd.read_parquet(cachedir / "predictions.parquet")
pfas_predictions = pfas_predictions.merge(chemicals[['name','chain_length','sulfonate']],on='name',how='left')
pfas_predictions['classification'] = pfas_predictions['chain_length'].apply(lambda x: "pfas_6+" if x >= 6 else "pfas <6")

# we want to take the classified neurotoxic chemicals and create a heatmap with them and the pfas
feature_path = "cache/notebooks/pfas_analysis/hepato-features.txt"
features = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines())

# region FIRST HEATMAP ======================================================
pdf = pfas_predictions[pfas_predictions['property_title'].isin(features)]
build_heatmap._generate_heatmap(
    pdf=pdf,
    output_path=cachedir / "heatmap.png",
    linecolor='black',
    fontcolor='white',
    dpi=600
)
# endregion

# region SHOW CHEMICALS WITH SHORT CHAINS THAT HAVE HIGH MEAN ACTIVITY VERSUS LOW
df = pfas_predictions[pfas_predictions['classification'] == 'pfas <6']
df = df.groupby(['name','inchi','chain_length','sulfonate'])['value'].mean().reset_index()
df = df.sort_values('value',ascending=False).reset_index(drop=True)

top3 = df.sort_values('value',ascending=False).iloc[[0,2,3,4]]['inchi']
bot3 = df.sort_values('value',ascending=True).iloc[0:4]['inchi']
top3_mols = [Chem.MolFromInchi(inchi) for inchi in top3]
bot3_mols = [Chem.MolFromInchi(inchi) for inchi in bot3]

# Combine into a grid image
all_mols = top3_mols + bot3_mols
legends = [f"Top {i+1}" for i in range(4)] + [f"Bot {i-2}" for i in range(4, 8)]
img = Draw.MolsToGridImage(all_mols,molsPerRow=8,subImgSize=(600, 300),legends=legends,useSVG=False)
img.save(str(cachedir / "topbot3.png"))
print(f"Saved top/bottom 5 molecule grid to {img_path}")

# endregion

# region MODEL DISTILLATION ======================================================
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.compose import ColumnTransformer

trainpred = pd.read_parquet("cache/projects/hepatotoxic/predictions.parquet")
trainpred = trainpred[trainpred['property_title'].isin(features)]
trainpred['classification'] = trainpred['classification'].apply(lambda x: 0 if x == 'Non-toxic' else 1)

# One-hot encode property_title, pivot to wide format for ML input
X = trainpred.pivot_table(index='name', columns='property_title', values='value', aggfunc='mean').fillna(0)
y = trainpred.drop_duplicates(subset='name').set_index('name').loc[X.index, 'classification']

# Train logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X, y)

# Predict
y_pred = model.predict_proba(X)[:, 1]
sort_ypred = pd.Series(y_pred, index=X.index, name='classification').sort_values(ascending=False)

# Save feature importances
feature_importance = pd.Series(model.coef_[0], index=X.columns)
important_features = feature_importance.abs().sort_values(ascending=False).head(20)

# Join predictions with metadata
pdf2 = pdf.copy()
pdf2 = pdf2.drop(columns=['classification'])
pdf2 = pdf2[pdf2['name'].isin(sort_ypred.index)]  # optional filter for clarity
pdf2 = pdf2.merge(sort_ypred, on='name', how='inner')

# Generate heatmap
build_heatmap._generate_heatmap_float(
    pdf=pdf2,
    output_path=cachedir / "heatmap.png",
    linecolor='black',
    fontcolor='white',
    dpi=600
)