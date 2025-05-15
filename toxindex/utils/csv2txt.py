import pathlib
import pandas as pd

cachedir = pathlib.Path('cache')
cachedir.mkdir(exist_ok=True)

project_name = 'nephrotoxic'
project_dir = cachedir / 'projects' / project_name

csv = pd.read_csv(project_dir / 'chemicals.csv', header=None)

txtfile = project_dir / 'chemicals.txt'
txtfile.write_text('\n'.join(csv.iloc[:, 0].astype(str)))
