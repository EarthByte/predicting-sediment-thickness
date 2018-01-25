# Predicting sediment thickness

Generate compacted sediment thickness and decompacted sediment rate grids for palaeo times using polynomials of ocean floor age and distance to passive margins.

## Workflow procedure

- Download the age grids (0-230Ma) and the topologies from ftp://ftp.earthbyte.org/Data_Collections/Muller_etal_2016_AREPS/
- Open the *generate_distance_grids.py* script and:
  + Set the *age_grid_dir* variable to the location of the downloaded age grids *netCDF_0-230Ma* directory.
  + Set the *topology_dir* variable to the location of the downloaded topologies *Muller_etal_AREPS_Supplement* directory.
  + Create the sub-directory *distances_1d*.
  + Set the *grid_spacing* variable to 1 (degree):
    * This takes about 8 hours on a 6-core (12-thread) system.
- Run the Python script:
    `python generate_distance_grids.py`
- The script outputs mean distance grids:
  + Have units in metres.
  + Are located in *distances_1d*.
- Open the *generate_predicted_sedimentation_grids.py* script and:
  + Set the *age_grid_dir* variable to the location of the downloaded age grids *netCDF_0-230Ma* directory.
  + Set the *distance_grid_dir* variable to the location of the generated distance grids (1 degree):
    * Currently this is *distances_1d*.
  + Create sub-directories *sedimentation_output/predicted_rate* and *sedimentation_output/predicted_thickness*.
  + Specify the grid spacing and time range variables:
    * Using 32-bit Python on Windows, the lowest grid spacing is 0.2 degrees (going lower will run out of memory):
      ^ This happens when using 32-bit pyGPlates - see http://www.gplates.org/docs/pygplates/pygplates_getting_started.html#using-the-correct-python-version
    * On 64-bit Mac and Linux you can go lower since pyGPlates (and hence Python) is 64-bit.
- Run the Python script:
    `python generate_predicted_sedimentation_grids.py`
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


## Miscellaneous

Global sediment thickness map:

Version 2...

https://www.ngdc.noaa.gov/mgg/sedthick/

Version 3...

One of the links at the bottom of this page...
https://doi.pangaea.de/10.1594/PANGAEA.835589?format=html#download

...seems to be a version of the global sediment thickness map, with their new grid merged into it
http://store.pangaea.de/Publications/WobbeF_et_al_2014/sedthick_world_v3_5min_epsg4326.nc 

