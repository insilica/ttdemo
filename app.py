import rapidfuzz
import numpy as np
import pathlib
import json

cachedir = pathlib.Path('cache')
cachedir.mkdir(exist_ok=True)

with open(cachedir / 'example_prediction.json', 'r') as f:
    propjson = json.loads(f.read())

proptitles = [prop.get('property').get('title') for prop in propjson]
propfile = cachedir / 'predicted_property_names.txt'
propfile.write_text('\n'.join(proptitles))

#%% GENERATE PROMPTS FOR TITLES
#
# ## cosmetics
# prompt_cosmetics = f"""I have selected chemicals that are actively used in the cosmetics industry.
# please select properties from the below list that would be interesting to toxicologists
# evaluating these compounds."""
# prompt_cosmetics = prompt_cosmetics + '\n\n' + '\n'.join(proptitles)
# prompt_cosmetics = prompt_cosmetics + '\n\n' + "pick ~100 of the most cosmetics relevant properties and output one per line and nothing else."
# prompt_cosmetics_file = cachedir / 'projects' / 'cosmetics' / 'property_prompt.txt'
# prompt_cosmetics_file.write_text(prompt_cosmetics)

#%% Fuzzy match relevant properties for each project to the predicted property names
def fuzzy_match_properties(input_path,output_path):
    relevant_properties = (input_path).read_text().splitlines()
    predicted_properties = (cachedir / 'predicted_property_names.txt').read_text().splitlines()
    #question = where did this come from? predicted_properties

    matched_properties = []
    # UV stability testing can be matched with chemical inhibition assay targeting timp2 which is incorrect
    for rel_prop in relevant_properties:
        # Calculate similarity scores between this relevant property and all predicted properties
        scores = [rapidfuzz.fuzz.token_sort_ratio(rel_prop.lower(), pred_prop.lower()) for pred_prop in predicted_properties]
        best_match_idx = int(np.argmax(scores))
        # matched_property = predicted_properties[best_match_idx].strip("[]")
        matched_properties.append(predicted_properties[best_match_idx])

    # write one per line to the project_dir / 'matched_properties.txt'
    output_path.write_text('\n'.join(matched_properties))
    return list(set(matched_properties))

projects = ['nephrotoxic']
# parse the chemicals

for project in projects:
    fuzzy_match_properties(
        input_path=cachedir / 'projects' / project / 'claude_relevant_properties.txt',
        output_path=cachedir / 'projects' / project / 'matched_properties.txt')

import toxindex.parse_chemicals as parse_chemicals
for project in projects:
    parse_chemicals.parse_chemicals(
        input_path=cachedir / 'projects' / project / 'chemicals.txt',
        output_path=cachedir / 'projects' / project / 'parsed_chemicals.csv'
    )

# categorize chemicals
import toxindex.categorize_chemicals as categorize_chemicals
for project in projects:
    categorize_chemicals.categorize_chemicals(
        input_path=cachedir / 'projects' / project / 'parsed_chemicals.csv',
        output_path=cachedir / 'projects' / project 
    )

# run predictions
import toxindex.predict_chemicals as predict_chemicals
for project in projects:
    predict_chemicals.predict_chemicals(
        input_path=cachedir / 'projects' / project / 'classified_chemicals.csv',
        output_path=cachedir / 'projects' / project / 'predictions.parquet'
    )

# build heatmaps
import toxindex.build_heatmap as build_heatmap
for project in projects:
    outdir = cachedir / 'projects' / project / 'heatmap_dir'
    outdir.mkdir(exist_ok=True)
    build_heatmap.build_heatmap(
        input_path=cachedir / 'projects' / project / 'predictions.parquet',
        output_path=outdir / 'heatmap.png'
    )

import toxindex.build_stripchart as build_stripchart
for project in projects:
    outdir = cachedir / 'projects' / project / 'stripchart_dir'
    outdir.mkdir(exist_ok=True)
    build_stripchart.build_stripchart(
        input_path=cachedir / 'projects' / project / 'predictions.parquet',
        output_path=outdir / 'stripchart.png'
    )
