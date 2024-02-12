
from ptt.utils.call_system_command import call_system_command
import os, shutil

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
""" ---------- Part 1 of the predicting-sediment-thickness workflow -----------
 This script creates grids (netcdfs) of the mean distance to passive margins through time

Requirements & Inputs:
    - Python  
    - Python scripts: ocean_basin_proximity.py, 
    - GMT 5 (or later)
    - Files associated with a tectonic model: agegrids, rotation file, plate boundaries, topologies

Outputs:
    - directory named 'distances_1d', with mean distance grids (in metres) through time.

2020-02-14: Added comments, removed hardcoded names (for rotation file, etc) from definitions
2022-08-29: Added more comments (NW)
"""

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------------------------------------
# --- Set paths and various parameters
# ------------------------------------------

# Use all CPUs.
#
# If False then use a single CPU.
# If True then use all CPUs (cores).
# If a positive integer then use that specific number of CPUs (cores).
#
#use_all_cpus = False
#use_all_cpus = 16
use_all_cpus = True

output_base_dir = '.'

# ------------------------------------------
# --- set times and spacing

# Generate distance grids for times in the range [min_time, max_time] at 'time_step' intervals.
min_time = 0
max_time = 250
time_step = 1

grid_spacing = 1

# The reference frame to generate the distance grids in.
#
# NOTE: The age grids must also be in this reference frame.
#
# Note: If the proximity/continent features have already been reconstructed
#       (eg, the features are actually time-dependent snapshots of reconstructions
#       generated from the continent contouring workflow) then they should remain
#       in the default reference frame (anchor plate zero in that workflow).
#       This is because they are assigned a plate ID of zero (in that workflow) and so
#       reconstructing them (in this workflow) relative to our anchor plate ID will
#       automatically position them correctly in our reference frame.
anchor_plate_id = 0

# ------------------------------------------
# --- input files
proximity_features_files = [
	'input_data/Global_EarthByte_GeeK07_COBLineSegments_2016_v4.gpmlz', # this is included in this repository
]

# Optional continent obstacles that the shortest distance path must go around (ie, water flowing around continents, rather than through).
# If not specifed then distances are minimum straight-line (great circle arc) distances from ocean points to proximity geometries.
# Obstacles can be both polygons and polylines.
#
#continent_obstacle_files = None
continent_obstacle_files = [
    '/Applications/GPlates_2.3.0/GeoData/FeatureCollections/Coastlines/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz',
]

# --- location of gplates files on your computer
# --- agegrids. Can be found here: https://www.earthbyte.org/gplates-2-3-software-and-data-sets/
#
#     The format string to generate age grid filenames (using the age grid paleo times in the range [min_time, max_time]).
#     Use a string section like "{:.1f}" to for the age grid paleo time. The ".1f" part means use the paleo time to one decimal place
#     (see Python\'s str.format() function) such that a time of 100 would be substituted as "100.0".
#     This string section will get replaced with each age grid time in turn (to generate the actual age grid filenames).
age_grid_filenames_format = '/Users/nickywright/Data/Age/Muller2019-Young2019-Cao2020_Agegrids/Muller2019-Young2019-Cao2020_netCDF/Muller2019-Young2019-Cao2020_AgeGrid-{:.0f}.nc'

# --- topologies and other files
data_dir = '/Applications/GPlates_2.3.0/GeoData/FeatureCollections/'
rotation_filenames = [
    '{}/Rotations/Muller2019-Young2019-Cao2020_CombinedRotations.rot'.format(data_dir)]

topology_filenames = [
    '{}/DynamicPolygons/Muller2019-Young2019-Cao2020_PlateBoundaries.gpmlz'.format(data_dir),
    '{}/DeformingLithosphere/Muller2019-Young2019-Cao2020_ActiveDeformation.gpmlz'.format(data_dir)]

# For each distance grid do not reconstruct ocean points earlier than 'max_topological_reconstruction_time'
# (each ocean point is reconstructed back to its age grid value or this value, whichever is smaller).
# This limit can be set to the earliest (max) reconstruction time of the topological model.
# If it's 'None' then only the age grid limits how far back each point is reconstructed.
max_topological_reconstruction_time = 1000  # can be None to just use age grid as the limit

# Optionally clamp mean distances to this value (in kms).
#clamp_mean_proximity_kms = None
clamp_mean_proximity_kms = 3000

output_dir = '{}/distances_{}d'.format(output_base_dir, grid_spacing)

# ------------------------------------------
# END USER INPUT
# ------------------------------------------

# -----
# make needed directories
if not os.path.exists(output_base_dir):
    os.makedirs(output_base_dir)


if not os.path.exists(output_dir):
    print('{} does not exist, creating now... '.format(output_dir))
    os.mkdir(output_dir)

# ----- 
def generate_distance_grids(times):
    py_cmd='python3'
    if os.environ.get('CONDA_PREFIX') or shutil.which('python3') is None:
        py_cmd = 'python'
    
    # Calling the ocean basin proximity script.
    command_line = [py_cmd, 'ocean_basin_proximity.py']

    # Rotation files.
    command_line.append('-r')
    command_line.extend('{}'.format(rotation_filename) for rotation_filename in rotation_filenames)

    # Proximity files.
    command_line.append('-m')
    command_line.extend('{}'.format(proximity_features_file) for proximity_features_file in proximity_features_files)
    # Proximity features are non-topological.
    command_line.append('-n')

    # If using continent obstacles.
    if continent_obstacle_files:
        command_line.append('--continent_obstacle_filenames')
        command_line.extend('{}'.format(continent_obstacle_file) for continent_obstacle_file in continent_obstacle_files)
    
    # Topological files.
    command_line.append('-s')
    command_line.extend('{}'.format(topology_filename) for topology_filename in topology_filenames)

    # Anchor plate ID.
    command_line.extend(['-a', '{}'.format(anchor_plate_id)])

    # Age grid filenames format.
    command_line.extend(['--age_grid_filenames_format', age_grid_filenames_format])
    # Age grid paleo times.
    command_line.append('--age_grid_paleo_times')
    command_line.extend(['{}'.format(time) for time in times])

    # Use all feature types in proximity file (according to Dietmar)...
    #command_line.extend(['-b', 'PassiveContinentalBoundary'])

    # Time increment is 1 Myr (this is for topological reconstruction of the ocean points).
    command_line.extend(['--time_increment', '1'])
    
    # If limiting the max topological reconstruction time.
    if max_topological_reconstruction_time is not None:
        command_line.extend(['-x', '{}'.format(max_topological_reconstruction_time)])
    
    # Grid spacing.
    command_line.extend(['-i', '{}'.format(grid_spacing)])

    # Optionally clamp mean proximity.
    if clamp_mean_proximity_kms:
        command_line.extend(['--clamp_mean_distance', str(clamp_mean_proximity_kms)])

    # Don't output distance grids for all reconstruction times.
    # Only outputting a single "mean" (over all reconstruction times) distance grid.
    #command_line.append('-d')

    # Output a "mean" (over all reconstruction times) distance grid.
    command_line.append('-j')

    # Generate a grd (".nc") file for each xyz file.
    command_line.append('-w')

    if use_all_cpus:
        # If 'use_all_cpus' is a bool (and therefore must be True) then use all available CPUs...
        if isinstance(use_all_cpus, bool):
            num_cpus = None  # use default of all available CPUs
        # else 'use_all_cpus' is a positive integer specifying the number of CPUs to use...
        elif isinstance(use_all_cpus, int) and use_all_cpus > 0:
            num_cpus = use_all_cpus
        else:
            raise TypeError('use_all_cpus: {} is neither a bool nor a positive integer'.format(use_all_cpus))
    else:
        num_cpus = 1
    
    # Number of cores.
    # If None then not specified, and defaults to using all available cores.
    if num_cpus:
        command_line.extend(['-c', '{}'.format(num_cpus)])

    # Distance grids output filename prefix.
    command_line.append('{}/distance_{}'.format(output_dir, grid_spacing))
    
    #print(' '.join(command_line))
    call_system_command(command_line)
    

    #
    # Rename the mean distance grids ('.nc') so that they start with 'mean_distance', and so that'time' is at the end
    # of the base filename (this way we can import them as time-dependent raster into GPlates version 2.0 and earlier).
    #
    # Also remove mean distance '.xy' files.
    #
    
    for time in  times:
        src_mean_distance_basename = '{}/distance_{}_{}_mean_distance'.format(output_dir, grid_spacing, time)
        dst_mean_distance_basename = '{}/mean_distance_{}d_{}'.format(output_dir, grid_spacing, time)
        
        # Rename '.nc' files.
        src_mean_distance_grid = src_mean_distance_basename + '.nc'
        dst_mean_distance_grid = dst_mean_distance_basename + '.nc'
        if os.access(dst_mean_distance_grid, os.R_OK):
            os.remove(dst_mean_distance_grid)
        if os.path.exists(src_mean_distance_grid):
            os.rename(src_mean_distance_grid, dst_mean_distance_grid)
        
        # Remove '.xy' files.
        src_mean_distance_xy = src_mean_distance_basename + '.xy'
        if os.access(src_mean_distance_xy, os.R_OK):
            os.remove(src_mean_distance_xy)


if __name__ == '__main__':

    #import time as time_prof

    print('Generating distance grids...')
    
    #tprof_start = time_prof.perf_counter()

    times = range(min_time, max_time + 1, time_step)
    #times = range(max_time, min_time - 1, -time_step) # Go backwards (can see results sooner).

    # Generate the distance grids.
    try:
        generate_distance_grids(times)
    except KeyboardInterrupt:
        pass
    
    tprof_end = time_prof.perf_counter()
    print(f"Total time: {tprof_end - tprof_start:.2f} seconds")
