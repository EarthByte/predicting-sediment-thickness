# README for updating

1. Set your agegrid and sediment thickness grid in `calculate_decompacted_sed_rate_and_xyz.sh`, and run.
** NOTE**: Make sure the script it calls (`average_sedimentation_rate.py`) is calculating (decompacted) rate.
This is in the `average_sedimentation_rate` function (around lines 273/274)

2. Run `update_alldata.ipynb`