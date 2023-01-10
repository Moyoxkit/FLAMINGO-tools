import pandas as pd
import swiftemulator as se
import numpy as np
import matplotlib.pyplot as plt
from velociraptor.observations import load_observations
from velociraptor.observations.objects import ObservationalData
import matplotlib as mpl
import emcee
import corner
import unyt
from unyt import Msun
from astropy.cosmology import WMAP9
from tqdm import tqdm
from sklearn.decomposition import IncrementalPCA
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
mpl.rcParams['figure.dpi'] = 200
import pandas as pd
import healpy as hp
from velociraptor import load as load_catalogue
import h5py as h5
from astropy.io import fits
from p_tqdm import p_map
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-o","--output", type=str)
parser.add_argument("-i","--input", type=str)
args = parser.parse_args()

big_cat = {}
for shell_nr in tqdm(range(4)):
    if shell_nr < 10:
        shell_name = "0" + str(shell_nr)
    else:
        shell_name = str(shell_nr)

    cat = pd.read_csv(args.input + "_" + shell_name + ".csv")

    if shell_nr ==0:
        for name in cat.keys():
            big_cat[name] = cat[name].values
    else:
        for name in cat.keys():
            all_prev = big_cat[name]
            new = cat[name].values
            added = np.append(all_prev,new)
            big_cat[name] = added

print(len(added))

out_cat = pd.DataFrame(big_cat)
out_cat.to_csv(args.output + ".csv")
