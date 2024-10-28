# Predicting sediment thickness

This workflow generates **compacted sediment thickness** and **decompacted sediment rate** grids for palaeo-times using polynomials of ocean floor age and distance to passive margins. The polynomial relationship is calibrated to present day and can be updated with new datasets if needed (see [below](#calculating-the-relationships-for-sedimentation-rate-and-thickness) for more information).

To generate sediment thickness and rate grids through time, all that is required is a GPlates-compatible plate motion model (specifically: rotation file(s), topology (or dynamic polgyon) file(s), and passive continental margin locations (COBs), and corresponding paleo-age grids. The latest plate model and time-evolving seafloor age grids can be downloaded [here](https://www.earthbyte.org/gplates-2-4-software-and-data-sets/).

> __Note:__ Currently you may notice artefacts in the mean-distance and sedimentation grids, such as alternating stripes, that are caused by non-optimal reconstruction of ocean basin points using the topological model. We plan to improve this in the future (by improving the collision detection of ocean basin points with topological plates and networks in pyGPlates).

## Dependencies

You'll also need to install the following Python dependencies:
* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [Generic Mapping Tools (GMT) ](https://www.generic-mapping-tools.org/)
* [PlateTectonicTools](https://github.com/EarthByte/PlateTectonicTools), or [GPlately](https://github.com/GPlates/gplately) (which now contains PlateTectonicTools).
* [PyGPlates](https://www.gplates.org/) version 0.30 or above (also is a dependency of `PlateTectonicTools` and `GPlately`).
* And, on Windows platforms, optionally install [psutil](https://pypi.org/project/psutil/) so that this workflow can use CPU cores in the *background* (ie, below-normal priority).


You can install these with conda:

```
conda create -n <conda-environment> -c conda-forge numpy scipy gmt platetectonictools
conda activate <conda-environment>
```

...where `<conda-environment>` should be replaced with the name of your conda environment.


## Releases
### v2.0
This release significantly __reduces__ running time and memory usage.

Other changes include:
- Can generate costly distance grids faster using a lower internal resolution (eg, 1 degree).
  - Then upscale to a higher resolution (eg, 0.1 degrees).
- Can specify maximum memory usage (to avoid out-of-memory crashes).
- Input parameters easier to configure.
  - Eg, specifying age grid filename format (where time is in the filename and how many decimal places).
- Can specify a non-zero reference frame (eg, if age grids are in a non-zero mantle frame).

### v1.1
This release contains the sediment thickness workflow with an updated calibration for sediment thickness and rate using [GlobSed](https://ngdc.noaa.gov/mgg/sedthick/) sediment thickness [(Straume et al. 2019)](https://doi.org/10.1029/2018GC008115) and age grid from 'Muller-2019-Young2019-Cao2020' in the [GPlates 2.3 sample data](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/) [(Zahirovic et al. 2022)](https://doi.org/10.1002/gdj3.146)

### v1.0
This release contains the original sediment thickness workflow from [Dutkiewicz et al. (2017)](https://doi.org/10.1002/2017GC007258).
The relationship for sedimentation rate and thickness was based on the calibration of the age grid from [Müller et al. (2016)](https://doi.org/10.1146/annurev-earth-060115-012211), and present-day sediment thickness of [Divins (2003)](https://www.ngdc.noaa.gov/mgg/sedthick/sedthick.html), incorporating additions by [Whittaker et al. (2013)](https://doi.org/10.1002/ggge.20181) and [Wobbe et al. (2014)](https://doi.org/10.1016/j.gloplacha.2014.09.006) for the Southern Ocean. 


## Workflow procedure

- Download paleo-age grids and associated topologies.
    - The latest plate model and paleo-age grids can be downloaded from [here](https://www.earthbyte.org/gplates-2-4-software-and-data-sets/).
    - Alternatively, Müller et al. (2016; AREPS) can be downloaded from [here](https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/)
- Open the `01_generate_distance_grids.py` script and:
    + Set the `min_time`, `max_time` and `time_step` time range variables for the times to generate distance grids.
    + Set the `age_grid_filenames_format` variable to the location/filenames of the downloaded age grids.
      + Note: This format string includes a pattern (such as `{:.1f}`) that will be substituted with the age grid paleo times.
    + Set the `data_dir` variable to the location of all your topological files.
    + Set the `rotation_filenames` and `topology_filenames` variables to match those in the `data_dir`.
    + Set the `max_topological_reconstruction_time` variable to oldest age supported by the topological model (eg, 250, 410, or 1000).
    + Set the `anchor_plate_id` variable to the reference frame in which to generate the distance grids.
      + Note: The age grids must also be in this reference frame.
    + Set the `proximity_features_files` variable to the passive margin files (that distances are calculated relative to).
      + Note: These can be passive margins generated from *contoured* continents (see [here](https://github.com/EarthByte/continent-contouring)).
    + Set the `continent_obstacle_files` variable to the continent files (obstacles to water flow).
      + The shortest distance (from ocean points to passive margins) must go around the continents (ie, water flows around continents).
      + This can be `None` to just use the minimum straight-line distance to passive margins (ignoring continent obstacles).
      + Note: These can be *contoured* continents (see [here](https://github.com/EarthByte/continent-contouring)).
      + Set the `plate_boundary_obstacles` variable to those plate boundary feature types that also act as obstacles (to water flow).
        + This should typically be left as the default (mid-ocean ridges and subduction zones), but you can change this if desired.
        + Note: This parameter is ignored unless `continent_obstacle_files` is also specified.
    + Set the `grid_spacing` variable to the desired grid spacing (in degrees, e.g.  0.1) of the generated distance grids.
      + The output distance grids are upscaled from the grid spacing used internally for computations (`internal_grid_spacing`).
    + Set the `internal_grid_spacing` variable to grid spacing (in degrees) used for internal distance computations.
      + Note: This parameter significantly affects the time it takes to generate distance grids in this workflow.
    + Set the `use_all_cpus` variable to the number of CPU cores to use (eg, False, True or a specific number).
    + Set the `max_memory_usage_in_gb` variable to the amount of memory (in GB) to use.
      + For example, set it to the amount of physical RAM (or less if running other workflows simultaneously).
- Run the Python script:
      `python 01_generate_distance_grids.py`
    + The script outputs mean distance grids:
        + Have units in metres.
        + Are located in *distances_`<grid_spacing>`d* (unles you've changed your output folder with `output_dir`).
    
    
- Open the `02_generate_predicted_sedimentation_grids.py` script and:
    + Set the `min_time`, `max_time` and `time_step` time range variables for the times to generate sedimentation grids.
    + Set the `age_grid_filenames_format` variable to the location/filenames of the downloaded age grids.
      + Note: This format string includes a pattern (such as `{:.1f}`) that will be substituted with the age grid paleo times.
    + Set the `distance_grid_spacing` variable to equal the `grid_spacing` variable used to generate the distance grids in part 1.
    + Set the `grid_spacing` variable to your desired spacing in degrees (of the output sedimentation grids).
    + Set the `use_all_cpus` variable to the number of CPU cores to use (eg, False, True or a specific number).
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

> __Note:__ The paper forgot to mention that mid-ocean ridges and subduction zones were hardwired to act as obstacles to water flow
(in addition to the continents). This has now been made a configurable parameter (see `plate_boundary_obstacles` above) that
defaults to mid-ocean ridges and subduction zones.
