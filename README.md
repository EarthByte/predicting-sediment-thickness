# Predicting sediment thickness

This workflow generates **compacted sediment thickness** and **decompacted sediment rate** grids for palaeo-times using polynomials of ocean floor age and distance to passive margins. The polynomial relationship is calibrated to present day and can be updated with new datasets if needed (see [below](#calculating-the-relationships-for-sedimentation-rate-and-thickness) for more information).

To generate sediment thickness and rate grids through time, all that is required is a GPlates-compatible plate motion model (specifically: rotation file(s), topology (or dynamic polgyon) file(s), and passive continental margin locations (COBs), and corresponding paleo-age grids. The latest plate model and time-evolving seafloor age grids can be downloaded [here](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/).

## Dependencies

You'll also need to install the following Python dependencies:
* [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools)

## Releases
### v1.1
This release contains the sediment thickness workflow with an updated calibration for sediment thickness and rate using [GlobSed](https://ngdc.noaa.gov/mgg/sedthick/) sediment thickness [(Straume et al. 2019)](https://doi.org/10.1029/2018GC008115) and age grid from 'Muller-2019-Young2019-Cao2020' in the [GPlates 2.3 sample data](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/) [(Zahirovic et al. 2022)](https://doi.org/10.1002/gdj3.146)

### v1.0
This release contains the original sediment thickness workflow from [Dutkiewicz et al. (2017)](https://doi.org/10.1002/2017GC007258).
The relationship for sedimentation rate and thickness was based on the calibration of the age grid from [Müller et al. (2016)](https://doi.org/10.1146/annurev-earth-060115-012211), and present-day sediment thickness of [Divins (2003)](https://www.ngdc.noaa.gov/mgg/sedthick/sedthick.html), incorporating additions by [Whittaker et al. (2013)](https://doi.org/10.1002/ggge.20181) and [Wobbe et al. (2014)](https://doi.org/10.1016/j.gloplacha.2014.09.006) for the Southern Ocean. 


## Workflow procedure

- Download paleo-age grids and associated topologies.
    - The latest plate model and paleo-age grids can be downloaded from [here](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/).
    - Alternatively, Müller et al. (2016; AREPS) can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/
- Open the `01_generate_distance_grids.py` script and:
    + Set the `data_dir` variable to the location of all your files (agegrids, topologies, etc)
    + Set the `agegrid_dir` variable to the location of the downloaded age grids.
    + Set the `agegrid_filename` and `agegrid_filename_ext` to correspond to the filename of the agegrids themselves, and the extension.
    + Set the `topology_dir` variable to the location of the downloaded topologies.
    + Set the `rotation_filenames` and `topology_filenames` to match those in the topology_dir.
    + Set the `grid_spacing` variable to your desired spacing in degrees, e.g.  1:
         * This takes about 8 hours on a 6-core (12-thread) system.
    + Open the script *ocean_basin_proximity.py* and set the *data_dir* and *coastline_filename* to match those in the model you're using.
- Run the Python script:
      `python 01_generate_distance_grids.py`
    + The script outputs mean distance grids:
        + Have units in metres.
        + Are located in *distances_${grid_spacing}d*. (or what you've named your output folder)
    
    
- Open the `02_generate_predicted_sedimentation_grids.py` script and:
    + Set the `agegrid_dir`, `agegrid_filename`, and agegrid_filename_ext` variables to correspond to the paleo-age grids directory.
    + Set the `distance_grid_dir` variable to the location of the generated distance grids from part 1.
    + Specify the `grid_spacing` and time range variables.
- Run the Python script:
    `python 02_generate_predicted_sedimentation_grids.py`
    + The script outputs predicted decompacted sedimentation rate grids:
        + Have units *cm/Ky* (not *m/My*).
        + Are located in *sedimentation_output/predicted_rate*.
    - The script also outputs predicted compacted sediment thickness grids:
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
