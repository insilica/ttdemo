import os
import biobricks as bb
import pandas as pd
import toxindex.utils.chemprop as chemprop
import pathlib

propjson = {}
with open('./prediction.json', 'r') as f:
    propjson = json.loads(f.read())

proptitles = df['title'].tolist()
# Create directory if it doesn't exist
os.makedirs('cache/projects', exist_ok=True)

# Write property titles to file, one per line
with open('cache/projects/property_names.txt', 'w') as f:
    for title in proptitles:
        f.write(f"{title}\n")



# 2. get grouped properties

# 3. build property lists for each project

# while I wait I will write a job posting for business development for insilica

benzene_inchi = 'InChI=1S/C6H6/c1-2-4-6-3-5-1/h1-6H'
prediction = chemprop.chemprop_predict_all(benzene_inchi)