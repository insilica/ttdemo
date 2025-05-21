import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif, RFE
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

import pathlib
import logging

logger = logging.getLogger(__name__)

def select_feature(input_path, output_path, method="random_forest", max_features=90):
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
    df = pd.read_parquet(input_path)

    # Load LLM suggested features
    feature_path = input_path.parent / 'matched_properties.txt'
    feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    feature_names = sorted(list(feature_names))
    df['is_in_lookup'] = df['property_title'].isin(feature_names)
    df = df[df['is_in_lookup']]

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

    # drop duplicate row
    class_series = pdf.drop_duplicates('name').set_index('name')['classification']
    # add classificaiton
    norm_values['substance_class'] = class_series

    # norm_values.to_csv('matrix.csv')
    # Extract features and labels
    X = norm_values.drop(columns=["substance_class"])
    y = LabelEncoder().fit_transform(norm_values["substance_class"])

    # method = "lasso" #15
    # method = "random_forest" #300
    # method = "mutual_info" #300
    # method = "rfe" #300

    if method == "lasso":
        model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)
        model.fit(X, y)
        selected_idx = np.where(model.coef_[0] != 0)[0]

    elif method == "random_forest":
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X, y)
        importances = model.feature_importances_
        selected_idx = np.argsort(importances)[::-1][:max_features]

    elif method == "mutual_info":
        mi = mutual_info_classif(X, y, random_state=42)
        selected_idx = np.argsort(mi)[::-1][:max_features]

    elif method == "rfe":
        base_model = LogisticRegression(solver='liblinear', max_iter=1000)
        rfe = RFE(estimator=base_model, n_features_to_select=max_features, step=0.1)
        rfe.fit(X, y)
        selected_idx = np.where(rfe.support_)[0]

    else:
        raise ValueError("Invalid method. Choose from: 'lasso', 'random_forest', 'mutual_info', 'rfe'")

    # Get column names for selected features
    selected_tokens = X.columns[selected_idx].tolist()
    selected_features = df[df['property_token'].isin(selected_tokens)]['property_title'].unique().tolist()
    print(f"[{method.upper():<15}] Selected {len(selected_idx):>3} features out of {len(feature_names)}")

    output_path.write_text('\n'.join(selected_features))
    return list(set(selected_features))

    




if __name__ == "__main__":
    project = 'hepatotoxic'

    cachedir = pathlib.Path('cache')
    cachedir.mkdir(exist_ok=True)

    outdir = cachedir / 'projects' / project 
    outdir.mkdir(exist_ok=True)

    input_path=cachedir / 'projects' / project / 'predictions.parquet'
    output_path=outdir / 'selected_properties.txt'
    select_feature(input_path, output_path)