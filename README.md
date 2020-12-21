# Molecular data explorer

Python script to create a web app with Dash to visualize molecular data and molecular geometries

## Requirements
- Python 3.6/3.7 

Use pip to install requirements:

`pip install -r requirements.txt`


## Usage

The app is executed as:

`python mol_view.py data.csv xyz_structures_directory`

With data_test.csv and test_structures provided here it should be:

`python mol_view.py data_test.csv test_structures`

Then open the following addres: http://127.0.0.1:8060/ in your browser

The data.csv will be read by pandas as `data = pd.read_csv('data.csv')`. 
'xyz_structures_directory' should be the directory containing the xyz structures in different files.

The data should be a file with coma separated values, with the first row containing the names of the columns. 
One column must be named 'xyz_file_name', for the app the load molecular structure corresponding to each row. 

For example:

data.csv------------------------ <br />
col_name1, col_name2, col_name3, col_name4, xyz_file_name  <br />
0, 0, 0, 0, file1.xyz  <br />
0, 0, 0, 0, file2.xyz  <br />
0, 0, 0, 0, file3.xyz  <br />
... <br />
end--------------------------------

as in data_test.csv
