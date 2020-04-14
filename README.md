# Predicting sediment thickness

Generate compacted sediment thickness and decompacted sediment rate grids for palaeo-times using polynomials of ocean floor age and distance to passive margins.

## Workflow procedure

- Download the age grids (0-230Ma) and associated topologies, e.g. from ftp://ftp.earthbyte.org/Data_Collections/Muller_etal_2016_AREPS/
- Open the *01_generate_distance_grids.py* script and:
    + Set the *data_dir* variable to the location of all your files (agegrids, topologies, etc)
    + Set the *age_grid_dir* variable to the location of the downloaded age grids *netCDF_0-230Ma* directory. Replace directory name if needed.
    + Set the *age_grid_filename* and *age_grid_filename_ext* to correspond to the filename of the agegrids themselves, and the extension.
    + Set the *topology_dir* variable to the location of the downloaded topologies *Muller_etal_AREPS_Supplement* directory.
    + Set the *rotation_filename*, *plateboundaries_filename*, and *topologybuildingblocks_filename* to match those in the topology_dir.
    + Set the *grid_spacing* variable to 1 (degree):
         * This takes about 8 hours on a 6-core (12-thread) system.
    + Open the script *ocean_basin_proximity.py* and set the *data_dir* and *coastline_filename* to match those in the model you're using.
- Run the Python script:
      `python 01_generate_distance_grids.py`
- The script outputs mean distance grids:
    + Have units in metres.
    + Are located in *distances_1d*. (or what you've named your output folder)
- Open the *02_generate_predicted_sedimentation_grids.py* script and:
    + Set the *age_grid_dir* variable to the location of the downloaded age grids *netCDF_0-230Ma* directory.
    + Set the *distance_grid_dir* variable to the location of the generated distance grids (1 degree):
        * Currently this is *distances_1d*.
    + Specify the grid spacing and time range variables:
        * Using 32-bit Python on Windows, the lowest grid spacing is 0.2 degrees (going lower will run out of memory):
          ^ This happens when using 32-bit pyGPlates - see http://www.gplates.org/docs/pygplates/pygplates_getting_started.html#using-the-correct-python-version
        * On 64-bit Mac and Linux you can go lower since pyGPlates (and hence Python) is 64-bit.
- Run the Python script:
    `python 02_generate_predicted_sedimentation_grids.py`
- The script outputs predicted sedimentation rate grids:
    + Have units *cm/Ky* (not *m/My*).
    + Are located in *sedimentation_output/predicted_rate*.
    + Also in that folder are compacted thickness grids obtained from the predicted rate via:
        * `thickness = compact(rate * age)`
- The script also outputs predicted sediment thickness grids:
    + Have units in metres.
    + Are located in *sedimentation_output/predicted_thickness*.
    + Also in that folder are rate grids obtained from the predicted thickness via:
        * `rate = decompact(thickness) / age`

## Calculating the relationships for sedimentation rate and thickness
The scripts to calculate the sedimentation rate and thickness relationships are in the folder `python_notebooks_and_input_data_archive`.
If a new relationship needs to be derived (for example, to be consistent with a different present-day age grid), run the `sediment_rate.ipynb` and `sediment_thick.ipynb` notebooks (with a modified `alldata` file, and update the script `generate_predicted_sedimentation_grids.py` to be consistent with the polynomial coefficients obtained from the jupyter notebooks.

## Miscellaneous

Most recent global sediment thickness map:

https://www.ngdc.noaa.gov/mgg/sedthick/

Previous versions: 

https://www.ngdc.noaa.gov/mgg/sedthick/sedthick.html

and 

One of the links at the bottom of this page...
https://doi.pangaea.de/10.1594/PANGAEA.835589?format=html#download

...seems to be a version of the global sediment thickness map, with their new grid merged into it
http://store.pangaea.de/Publications/WobbeF_et_al_2014/sedthick_world_v3_5min_epsg4326.nc 


## Reference

Dutkiewicz, A., Müller, R.D., Wang, X., O’Callaghan, S., Cannon, J. and Wright, N.M., 2017. Predicting sediment thickness on vanished ocean crust since 200 Ma. Geochemistry, Geophysics, Geosystems, 18, 4586-4603, DOI:  https://doi.org/10.1002/2017GC007258.
