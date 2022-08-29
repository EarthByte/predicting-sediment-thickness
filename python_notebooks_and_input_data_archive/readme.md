# README for updating calibration files for the sediment thickness and rate relationships

To update the sediment thickness and rate relationships, you first need to update `alldata`. 

`alldata` is essentially a text file with the columns:
lon | lat | sedrate | sedthick | age | PCB_mean_dist | river_mouth: amazon | river_mouth: congo | river_mouth: ganges | river_mouth: godavari | river_mouth: Indus | river_mouth: magdalena | river_mouth: mahanadi | river_mouth: mississippi | river_mouth: narmada | river_mouth: niger | river_mouth: orinoco | river_mouth: paleo-congo | river_mouth: parana

This can be done by:

1. Seting your preferred age grid and sediment thickness grid in `calculate_decompacted_sed_rate_and_xyz.sh`, and run. 
    * **NOTE**: Make sure the script it calls (`average_sedimentation_rate.py`) is calculating (decompacted) rate.
This is in the `average_sedimentation_rate` function (around lines 273/274). 
    * This will also create 1 degree versions of your age and sediment thickness grids.
2. Launch the notebook `update_alldata.ipynb`. This notebook:
    * Copies your original `alldata` file to `alldata_orig`
    * Reads in the new 1° versions of the age, sediment thickness, and sediment rate grids (created in part 1).
    * Replaces the age/sediment thickness/sediment rate columns of `alldata` with these new grids
    * Saves out the data as `alldata`
