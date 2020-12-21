# Molecular data explorer

Python script to create a web app with Dash to visualize molecular data and 3d molecular geometries


## Requirements
- Python 3.6/3.7 

Use pip to install requirements:

`pip install -r requirements.txt`


## Usage

The app is executed as:

`python mol_view.py data.csv xyz_structures_directory`

The data.csv will be read by pandas as `data = pd.read_csv('data.csv')`
The format should be:

data.csv------------------------
col_name1, col_name2, col_name3, ... <br />
val1, val2, val3,...
...
--------------------------------

