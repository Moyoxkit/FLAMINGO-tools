"""
This script takes the output from the halo lightcone (created using the black hole lightcone) to some select quantities from the SOAP catalogue.
As this requires the loading of the SOAP catalogue at each redshift this is potentially slow.

This script uses the following inputs

-o base name of the output file
-z maximal redshift used
-l folder which has all the lightcone files
-s folder that has all the soap files
"""

import pandas as pd
import numpy as np
import healpy as hp
import unyt
from unyt import Msun
from tqdm import tqdm
import h5py as h5
import argparse
from velociraptor import load as load_catalogue

parser = argparse.ArgumentParser()
parser.add_argument("-o","--output", type=str)
parser.add_argument("-z","--redshift", type=float)
parser.add_argument("-l","--lightcone",type=str)
parser.add_argument("-s","--soap",type=str)
args = parser.parse_args()

# Calculate the number of redshift slices
max_index = int(np.rint(args.redshift / 0.05))

for shell_nr in tqdm(range(max_index)):
    if shell_nr < 10:
        shell_name = "0" + str(shell_nr)
    else:
        shell_name = str(shell_nr)

    catalogue = h5.File(args.lightcone + "/lightcone_halos_00" + str(shell_name) + ".hdf5","r")

    # Load the correct attribute and extract the halo type
    SH_cat = catalogue["Subhalo"]
    halotype = SH_cat["Structuretype"][:]

    # Make sure we always use the correct snapshot
    Snapnum = SH_cat["SnapNum"][:][halotype==10]

    # For the 2.8 and 1Gpc 77 should be replaced by 78 as we have one more snapshot
    Snap_num = 77 - shell_nr

    #Make a mask for the Snapnum (This is mostly a sanity check as all halos in the lightcone have the same Snapnum)
    Snapmask = Snapnum == int(Snap_num)

    #Mask out the very small objects
    M_fof = SH_cat["Mass_tot"][:][halotype==10][Snapmask]*1e10
    M_fof_mask = M_fof > 10**(13.0)

    # Load the Halo indices used for the matching
    lc_IDs = SH_cat["ID"][:][halotype==10][Snapmask][M_fof_mask]

    # Load the coordinates and redshift, and calculate the RA and dec
    xcoords = SH_cat["Xcmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    ycoords= SH_cat["Ycmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    zcoords = SH_cat["Zcmbp_bh"][:][halotype==10][Snapmask][M_fof_mask]
    z_of_cluster = SH_cat["Redshift"][:][halotype==10][Snapmask][M_fof_mask]
    lc_coords = np.array([xcoords,ycoords,zcoords])
    theta_c, phi_c = hp.rotator.vec2dir(lc_coords,lonlat=True)
    M_fof = M_fof[M_fof_mask]

    # Load the outputs into a dictonary to be saved as a csv
    to_save = {}
    to_save["id"] = lc_IDs
    to_save["redshift"] = z_of_cluster
    to_save["theta_on_lc"] = theta_c
    to_save["phi_on_lc"] = phi_c
    to_save["M_fof_lc"] = M_fof
    to_save["x_lc"] = xcoords
    to_save["y_lc"] = ycoords
    to_save["z_lc"] = zcoords

    # Now we need the to link the SOAP catalog
    catalogue = h5.File(args.soap + "/halo_properties_00" + str(Snap_num) + ".hdf5","r")

    # Start by loading the needed quantities
    hosthaloid = catalogue["VR/HostHaloID"][:]
    normalid = catalogue["VR/ID"][:]
    # Having Mfof from SOAP and the halo lightcone is a good check to see if the matching worked
    SOAP_Mfof = catalogue["FOFSubhaloProperties/TotalMass"][:]

    # This part will depend on what you personally need so edit here
    SOAP_M500 = catalogue["SO/500_crit/TotalMass"][:]
    SOAP_CY = catalogue["SO/5xR_500_crit/ComptonY"][:]
    SOAP_Mfof = catalogue["FOFSubhaloProperties/TotalMass"][:]
    SOAP_XR0 = catalogue["SO/500_crit/XRayLuminosity"][:,0]
    SOAP_XR1 = catalogue["SO/500_crit/XRayLuminosity"][:,1]
    SOAP_XR2 = catalogue["SO/500_crit/XRayLuminosity"][:,2]

    # We take advance of the fact that for SOAP IDs and Indices are seperated by 1 to fetch the correct halos from SOAP
    one_line_get = np.array(lc_IDs - 1,dtype=int)

    # Take only the matched indices and store them in the dictionary
    to_save["M500"] = SOAP_M500[one_line_get]
    to_save["Compton_y"] = SOAP_CY[one_line_get]
    to_save["Xray band 0"] = SOAP_XR0[one_line_get]
    to_save["Xray band 1"] = SOAP_XR0[one_line_get]
    to_save["Xray band 2"] = SOAP_XR0[one_line_get]
    to_save["M_fof_soap"] = SOAP_Mfof[one_line_get]

    # Save the final dictionary to csv
    out_cat = pd.DataFrame(to_save)
    out_cat.to_csv(args.output + "_" + shell_name + ".csv")
