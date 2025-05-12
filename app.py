import os
import rapidfuzz
import numpy as np
import pandas as pd
import biobricks as bb
import pandas as pd
import toxindex.utils.chemprop as chemprop
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

## cosmetics
prompt_cosmetics = f"""I have selected chemicals that are actively used in the cosmetics industry.
please select properties from the below list that would be interesting to toxicologists 
evaluating these compounds."""
prompt_cosmetics = prompt_cosmetics + '\n\n' + '\n'.join(proptitles)
prompt_cosmetics = prompt_cosmetics + '\n\n' + "pick ~100 of the most cosmetics relevant properties and output one per line and nothing else."
prompt_cosmetics_file = cachedir / 'projects' / 'cosmetics' / 'property_prompt.txt'
prompt_cosmetics_file.write_text(prompt_cosmetics)

## food-coloring
prompt_food_coloring = f"""I have selected chemicals that are actively used for food-coloring. 
please select properties from the below list that would be interesting to toxicologists 
evaluating these compounds."""
prompt_food_coloring = prompt_food_coloring + '\n\n' + '\n'.join(proptitles)
prompt_food_coloring = prompt_food_coloring + '\n\n' + "pick ~100 of the most food-coloring relevant properties and output one per line and nothing else."

prompt_food_coloring_file = cachedir / 'projects' / 'food-coloring' / 'property_prompt.txt'
prompt_food_coloring_file.write_text(prompt_food_coloring)

## food-contact-materials
prompt_food_contact = f"""I have selected chemicals that are actively used for food-contact-materials. 
please select properties from the below list that would be interesting to toxicologists 
evaluating these compounds."""
prompt_food_contact = prompt_food_contact + '\n\n' + '\n'.join(proptitles)
prompt_food_contact = prompt_food_contact + '\n\n' + "pick ~100 of the most food-contact-materials relevant properties and output one per line and nothing else."
prompt_food_contact_file = cachedir / 'projects' / 'food-contact-materials' / 'property_prompt.txt'
prompt_food_contact_file.write_text(prompt_food_contact)

## pfas
prompt_pfas = f"""I have selected per- and polyfluoroalkyl substances (PFAS) chemicals. 
Please select properties from the below list that would be interesting to toxicologists 
evaluating these compounds."""
prompt_pfas = prompt_pfas + '\n\n' + '\n'.join(proptitles)
prompt_pfas = prompt_pfas + '\n\n' + "pick ~100 of the most PFAS relevant properties and output one per line and nothing else."
prompt_pfas_file = cachedir / 'projects' / 'pfas' / 'property_prompt.txt'
prompt_pfas_file.write_text(prompt_pfas)

#%% Fuzzy match relevant properties for each project to the predicted property names
def fuzzy_match_properties(project_name):
    relevant_properties = (cachedir / 'projects' / project_name / 'claude_relevant_properties.txt').read_text().splitlines()
    predicted_properties = (cachedir / 'predicted_property_names.txt').read_text().splitlines()
    
    matched_properties = []
    for rel_prop in relevant_properties:
        # Calculate similarity scores between this relevant property and all predicted properties
        scores = [fuzz.token_sort_ratio(rel_prop.lower(), pred_prop.lower()) for pred_prop in predicted_properties]
        best_match_idx = int(np.argmax(scores))
        matched_properties.append(predicted_properties[best_match_idx])
    
    # write one per line to the project_dir / 'matched_properties.txt'
    matched_properties_file = cachedir / 'projects' / project_name / 'matched_properties.txt'
    matched_properties_file.write_text('\n'.join(matched_properties))
    return list(set(matched_properties))

for project in ['cosmetics', 'food-coloring', 'food-contact-materials', 'pfas']:
    fuzzy_match_properties(project)

#%% Build a heatmap for each project from the chemicals.txt and matched_properties.txt
# 1. parse the chemicals to inchi
# 2. run the chemprop predict_all function on each chemical
# 3. filter to the matched properties
# 4. build a row and column clusterd heatmap of the chemicals and the properties
# 5. generate a ppt slide for each project

