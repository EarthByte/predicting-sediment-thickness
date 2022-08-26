#!/bin/bash
# script to run average_sedimentation_rate.py

sedfile=/Users/nickywright/Data/SedimentThickness/GlobSed_package2/GlobSed-v2.nc
agegrid=/Users/nickywright/Data/Age/Muller2019-Young2019-Cao2020_Agegrids/Muller2019-Young2019-Cao2020_netCDF/Muller2019-Young2019-Cao2020_AgeGrid-0.nc


python average_sedimentation_rate.py -m $sedfile -g $agegrid -i 0.1 -w -- average_sedimentation_rate

# created the following decompacted sedimentation thickness grid by editing line 271/272 of the average_sedimentation script
# e.g. one line divides by age to get rate, the other just gets decompacted sediment thickness
# python average_sedimentation_rate.py -m $sedfile -g $agegrid -i 0.1 -w -- decompacted_sediment_thickness


# convert age and sediment thickness grids to 1d
echo "... converting grids to 1degree"

gmt grdsample $sedfile -I1d -Gsedgrid_1d.nc -V
gmt grd2xyz sedgrid_1d.nc > sedgrid_1d.xyz

gmt grdsample $agegrid -I1d -Gagegrid_1d.nc -V
gmt grd2xyz agegrid_1d.nc > agegrid_1d.xyz

gmt grdsample average_sedimentation_rate.grd -I1d -Gaverage_sedimentation_rate_1d.nc -V