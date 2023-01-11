#!/bin/bash

# Run the linking script on the 1Gpc to z=0.6
python3 link_halo_lightcone_to_SOAP.py -o test -z 0.6 -l /cosma8/data/dp004/jch/FLAMINGO/ScienceRuns/L1000N1800/HYDRO_FIDUCIAL/lightcone_halos/lightcone0 -s /cosma8/data/dp004/flamingo/Runs/L1000N1800/HYDRO_FIDUCIAL/SOAP

# Combine the resulting csvs
python3 combine_linked_lcSOAP_catalogues.py -i test -o test -c
