import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif, RFE
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

import pathlib
import logging

logger = logging.getLogger(__name__)

def select_feature(input_path, output_path, max_features=1000):
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
    # feature_path = input_path.parent / 'matched_properties.txt'
    # feature_names = set(line.strip() for line in pathlib.Path(feature_path).read_text().splitlines() if line.strip())
    # feature_names = sorted(list(feature_names))
    feature_path = input_path.parent / 'chatgpt_selected_features.csv'
    feature_names = pd.read_csv(feature_path)
    # feature_names = feature_names[feature_names['keyword_occurrence'] > 0]
    feature_names = feature_names['property_title']


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
    methods = ["lasso","random_forest","mutual_info","rfe"]
    for method in methods:

        # lasso: Selects features with non-zero coefficients from L1-regularized logistic regression; 
        if method == "lasso":
            model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)
            model.fit(X, y)
            selected_idx = np.where(model.coef_[0] != 0)[0]
            scores = np.abs(model.coef_[0])[selected_idx] #score = absolute coefficient value

        # random_forest: Selects top features by tree-based impurity reduction; 
        elif method == "random_forest":
            model = RandomForestClassifier(n_estimators=100, random_state=42)
            model.fit(X, y)
            importances = model.feature_importances_
            selected_idx = np.argsort(importances)[::-1][:max_features]
            scores = importances[selected_idx] #score = feature importance from Random Forest

        # mutual_info: Selects features most informative about the target;
        elif method == "mutual_info":
            mi = mutual_info_classif(X, y, random_state=42)
            selected_idx = np.argsort(mi)[::-1][:max_features]
            scores = mi[selected_idx] #score = mutual information with the class label

        # rfe: Recursively removes features based on importance;
        elif method == "rfe":
            base_model = LogisticRegression(solver='liblinear', max_iter=1000)
            rfe = RFE(estimator=base_model, n_features_to_select=max_features, step=0.1)
            rfe.fit(X, y)
            selected_idx = np.where(rfe.support_)[0]
            #score = coefficient value from Logistic Regression
            scores = rfe.estimator_.coef_[0] 
            

        else:   
            raise ValueError("Invalid method. Choose from: 'lasso', 'random_forest', 'mutual_info', 'rfe'")

        scores = scores.round(4)


        # Step 1: Get selected tokens and scores
        selected_property_tokens = np.array(X.columns[selected_idx])
        scores = np.array(scores)

        # Step 2: Create a clean token-to-title mapping
        token_title_map = df.drop_duplicates(subset="property_token").set_index("property_token")["property_title"]

        # Step 3: Filter to tokens that exist in the map
        valid_mask = np.isin(selected_property_tokens, token_title_map.index)
        valid_tokens = selected_property_tokens[valid_mask]
        valid_scores = scores[valid_mask]
        valid_titles = token_title_map.loc[valid_tokens].values

        # Step 4: Assemble DataFrame
        temp_df = pd.DataFrame({
            "property_token": valid_tokens,
            "property_title": valid_titles,
            "score": valid_scores
        })

        # Step 5: Drop duplicates by title (optional: keep highest-scoring one)
        temp_df = temp_df.sort_values("score", ascending=False).drop_duplicates("property_title")
        temp_df = temp_df.sort_values("score", ascending=False).reset_index(drop=True)
        temp_df["rank"] = np.arange(1, len(temp_df) + 1)

        # Final result
        output_df = temp_df[["rank", "score", "property_title", "property_token"]]

        # Sort by score descending and reset rank
        output_df = output_df.sort_values(by="score", ascending=False).reset_index(drop=True)
        output_df["rank"] = np.arange(1, len(output_df) + 1)

        # compare the rank vs LLM keyword matching score
        # keyword_df = pd.read_csv(input_path.parent / 'chatgpt_selected_features.csv')
        # keyword_df = keyword_df.drop_duplicates(subset='property_title')
        # keyword_df = keyword_df[['property_title','keyword_occurrence']]
        # output_df = output_df.merge(keyword_df, on='property_title',how='left')
        # output_df = output_df.dropna(subset=['keyword_occurrence'])

        print(f"[{method.upper():<15}] Selected {len(output_df):>3} features out of {len(feature_names)}")


        output_csv_path = output_path / f"{method}_selected_properties.csv"
        output_df.to_csv(output_csv_path, index=False)

        # Also save plain list as TXT if needed
        output_txt_path = output_path / f"{method}_selected_properties.txt"
        output_txt_path.write_text('\n'.join(output_df.property_title.to_list()))

    return output_df


if __name__ == "__main__":
    project = 'hepatotoxic'

    cachedir = pathlib.Path('cache')
    cachedir.mkdir(exist_ok=True)

    outdir = cachedir / 'projects' / project / 'selected_properties'
    outdir.mkdir(exist_ok=True)

    input_path=cachedir / 'projects' / project / 'predictions.parquet'
    output_path=outdir 
    select_feature(input_path, output_path)

# # Step 1: Count how many unique tokens map to each title
# title_to_token_counts = df.groupby("property_title")["property_token"].nunique().reset_index()
# # Step 2: Filter to only those titles with more than one token
# duplicate_titles = title_to_token_counts[title_to_token_counts["property_token"] > 1]
