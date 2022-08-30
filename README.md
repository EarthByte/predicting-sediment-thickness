# Predicting sediment thickness

Generate **compacted sediment thickness** and **decompacted sediment rate** grids for palaeo-times using polynomials of ocean floor age and distance to passive margins.

## Releases
### v1.1
This release contains the sediment thickness workflow with an updated calibration for sediment thickness and rate using [GlobSed](https://ngdc.noaa.gov/mgg/sedthick/) sediment thickness [(Straume et al. 2019)](https://doi.org/10.1029/2018GC008115) and age grid from 'Muller-2019-Young2019-Cao2020' in the [GPlates 2.3 sample data](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/) [(Zahirovic et al. 2021)](https://doi.org/10.1002/gdj3.146)

### v1.0
This release contains the original sediment thickness workflow from [Dutkiewicz et al. (2017)](https://doi.org/10.1002/2017GC007258).
The relationship for sedimentation rate and thickness was based on the calibration of the age grid from [Müller et al. (2016)](https://doi.org/10.1146/annurev-earth-060115-012211), and present-day sediment thickness of [Divins (2003)](https://www.ngdc.noaa.gov/mgg/sedthick/sedthick.html), incorporating additions by [Whittaker et al. (2013)](https://doi.org/10.1002/ggge.20181) and [Wobbe et al. (2014)](https://doi.org/10.1016/j.gloplacha.2014.09.006) for the Southern Ocean. 


## Workflow procedure

- Download the age grids (0-230Ma) and associated topologies, e.g. from the webdav folder https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/ (Google how to connect to a webdav folder. It is nice and easy) or if you prefer the old FTP protocal ftp://ftp.earthbyte.org/Data_Collections/Muller_etal_2016_AREPS/ 
- Open the *01_generate_distance_grids.py* script and:
    + Set the *data_dir* variable to the location of all your files (agegrids, topologies, etc)
    + Set the *age_grid_dir* variable to the location of the downloaded age grids *netCDF_0-230Ma* directory. Replace directory name if needed.
    + Set the *age_grid_filename* and *age_grid_filename_ext* to correspond to the filename of the agegrids themselves, and the extension.
    + Set the *topology_dir* variable to the location of the downloaded topologies *Muller_etal_AREPS_Supplement* directory.
    + Set the *rotation_filenames* and *topology_filenames* to match those in the topology_dir.
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
- The script also outputs predicted sediment thickness grids:
    + Have units in metres.
    + Are located in *sedimentation_output/predicted_thickness*.

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
