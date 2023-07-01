# Title
Data Analysis Functions for Immuno-Fluorescence Data derived from Omero-Screen

## Status
In Development

## Authors
Helfrid Hochegger, Robert Zach, Ryan Yue

### Dependencies
Requires Python 3.6 or greater 

## Contact
Created by Helfrid Hochegger
email: hh65@sussex.ac.uk

## Licence
MIT Licence

## Usage

1. Count Module
To generate a simple count of the number of cells in each well, use the following command:
from ifanalysis.counts import count_cells
count_cells(df, conditions, 'CTR', path=OUTPUT_PATH)

df is the dataframe containing the data from the omero_Screen run
conditions is a list of the conditions to be analysed
'CTR' is the name of the control condition to compare the other conditions to
OUTPUT_PATH is the path to the folder where the output should be saved
example figure: 

2. Cell Cycle Module
To analyse cell cycle data, use the following command

from ifanalysis.cellcycle import standard_cellcycleplots
standard_cellcycleplots(df, conditions, file_path=OUTPUT_PATH))

df is the dataframe containing the data from the omero_Screen run
conditions is a list of the conditions to be analysed
OUTPUT_PATH is the path to the folder where the output should be saved