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
args = parser.parse_args()

returned_names = ["id", "M500", "Compton_y", "Xray band 0", "Xray band 1", "Xray band 2"]

for shell_nr in tqdm(range(4)):
    if shell_nr < 10:
        shell_name = "0" + str(shell_nr)
    else:
        shell_name = str(shell_nr)
    catalogue = h5.File("/cosma8/data/dp004/jch/FLAMINGO/ScienceRuns/L1000N1800/HYDRO_FIDUCIAL/lightcone_halos/lightcone0/lightcone_halos_00" + str(shell_name) + ".hdf5","r")
    
    SH_cat = catalogue["Subhalo"]
    halotype = SH_cat["Structuretype"][:]

    xcoords = SH_cat["Xcmbp_bh"][:][halotype==10]
    ycoords= SH_cat["Ycmbp_bh"][:][halotype==10]
    zcoords = SH_cat["Zcmbp_bh"][:][halotype==10]
    lc_coords = np.array([xcoords,ycoords,zcoords])

    #Apply both masks to the fields we want to laod
    Snapnum = SH_cat["SnapNum"][:][halotype==10]

    Snap_num = 77 - shell_nr
    #Make a mask for the Snapnum
    Snapmask = Snapnum == int(Snap_num)

    #Mask out the very small objects
    M_fof = SH_cat["Mass_tot"][:][halotype==10][Snapmask]*1e10
    M_fof_mask = M_fof > 10**(13.0)

    lc_IDs = SH_cat["ID"][:][halotype==10][Snapmask][M_fof_mask]
    xcoords = SH_cat["Xcmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    ycoords= SH_cat["Ycmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    zcoords = SH_cat["Zcmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    z_of_cluster = SH_cat["Redshift"][:][halotype==10][Snapmask][M_fof_mask]
    lc_coords = np.array([xcoords,ycoords,zcoords])
    theta_c, phi_c = hp.rotator.vec2dir(lc_coords,lonlat=True)
    M_fof = M_fof[M_fof_mask]

    catalogue = h5.File("/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP/halo_properties_00" + str(Snap_num) + ".hdf5","r")
    vrcar = load_catalogue("/cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/VR/catalogue_00" + str(Snap_num) + "/vr_catalogue_00" + str(Snap_num) +".properties.0","r")

    hosthaloid = np.array(vrcar.ids.hosthaloid.value,dtype=int)
    normalid = np.array(vrcar.ids.id.value,dtype=int)

    SOAP_M500 = catalogue["SO/500_crit/TotalMass"][:]
    SOAP_CY = catalogue["SO/5xR_500_crit/ComptonY"][:]
    SOAP_Mfof = catalogue["FOFSubhaloProperties/TotalMass"][:]
    SOAP_XR0 = catalogue["SO/500_crit/XRayLuminosity"][:,0]
    SOAP_XR1 = catalogue["SO/500_crit/XRayLuminosity"][:,1]
    SOAP_XR2 = catalogue["SO/500_crit/XRayLuminosity"][:,2]

    one_line_get = np.array(lc_IDs - 1,dtype=int)

    to_save = {}

    to_save["id"] = lc_IDs
    to_save["M500"] = SOAP_M500[one_line_get]
    to_save["Compton_y"] = SOAP_CY[one_line_get]
    to_save["Xray band 0"] = SOAP_XR0[one_line_get]
    to_save["Xray band 1"] = SOAP_XR0[one_line_get]
    to_save["Xray band 2"] = SOAP_XR0[one_line_get]
    to_save["redshift"] = z_of_cluster
    to_save["theta_on_lc"] = theta_c
    to_save["phi_on_lc"] = phi_c
    to_save["M_fof_lc"] = M_fof
    to_save["M_fof_soap"] = SOAP_Mfof[one_line_get]
    to_save["x_lc"] = xcoords
    to_save["y_lc"] = ycoords
    to_save["z_lc"] = zcoords


    out_cat = pd.DataFrame(to_save)
    out_cat.to_csv(args.output + "_" + shell_name + ".csv")
