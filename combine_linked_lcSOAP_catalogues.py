"""
Script to combine the outputs 

this has three inputs
-i the same name used as output for the linking script
-o the name of the output file (not including the .csv part)
-c A flag to delete the input files when the final catalog is saved
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
from glob import glob
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-o","--output", type=str)
parser.add_argument("-i","--input", type=str)
parser.add_argument("-c","--clean",action='store_true')

# Parse the folder for the individual files
args = parser.parse_args()
files = [Path(x) for x in glob(args.input + "*.csv")]

# Loop over the different redshifts and concatenate into a big dictonary
big_cat = {}
for shell_nr, file in tqdm(enumerate(files)):
    cat = pd.read_csv(file)
    if shell_nr == 0:
        for name in cat.keys():
            big_cat[name] = cat[name].values
    else:
        for name in cat.keys():
            all_prev = big_cat[name]
            new = cat[name].values
            added = np.append(all_prev,new)
            big_cat[name] = added

# Save the output
out_cat = pd.DataFrame(big_cat)
out_cat.to_csv(args.output + ".csv")

# Remove the individual files if requested
if args.clean:
    for file in files:
        os.remove(file)
