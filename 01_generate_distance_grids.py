
from ptt.utils.call_system_command import call_system_command
import math
import multiprocessing
import os, shutil
import sys

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

output_base_dir = '.'

# ------------------------------------------
# --- set times and spacing

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

# --- location of gplates files on your computer
# DON'T FORGET TO UPDATE ocean_basin_proximity.py!
# --- agegrids. Can be found here: https://www.earthbyte.org/gplates-2-3-software-and-data-sets/
agegrid_dir = '/Users/nickywright/Data/Age/Muller2019-Young2019-Cao2020_Agegrids/Muller2019-Young2019-Cao2020_netCDF'   # change folder name if needed
agegrid_filename_prefix = 'Muller2019-Young2019-Cao2020_AgeGrid-'    # everything before 'time'
agegrid_filename_suffix = ''    # everything after 'time' (excluding the extension), eg, "Ma"
agegrid_filename_ext = 'nc'   # generally 'nc', but sometimes is 'grd'. Do not include the period

# --- topologies and other files
data_dir = '/Applications/GPlates_2.3.0/GeoData/FeatureCollections/'
rotation_filenames = [
    '{}/Rotations/Muller2019-Young2019-Cao2020_CombinedRotations.rot'.format(data_dir)]

topology_filenames = [
    '{}/DynamicPolygons/Muller2019-Young2019-Cao2020_PlateBoundaries.gpmlz'.format(data_dir),
    '{}/DeformingLithosphere/Muller2019-Young2019-Cao2020_ActiveDeformation.gpmlz'.format(data_dir)]

proximity_threshold_kms = 3000

output_dir = '{}/distances_{}d'.format(output_base_dir, grid_spacing)

# Use all CPUs.
#
# If False then use a single CPU.
# If True then use all CPUs (cores).
# If a positive integer then use that specific number of CPUs (cores).
#
#use_all_cpus = False
#use_all_cpus = 4
use_all_cpus = True

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
def generate_distance_grid(time):
    py_cmd='python3'
    if os.environ.get('CONDA_PREFIX') or shutil.which('python3') is None:
        py_cmd = 'python'
    
    command_line = [py_cmd, 'ocean_basin_proximity.py']
    command_line.extend(['-r'])
    command_line.extend('{}'.format(rotation_filename) for rotation_filename in rotation_filenames)
    command_line.extend(['-m'])
    command_line.extend('{}'.format(proximity_features_file) for proximity_features_file in proximity_features_files)
    command_line.extend(['-s'])
    command_line.extend('{}'.format(topology_filename) for topology_filename in topology_filenames)
    command_line.extend(['-a'])
    command_line.extend(['{}'.format(anchor_plate_id)])
    command_line.extend([
            '-g',
            '{}/{}{}{}.{}'.format(agegrid_dir, agegrid_filename_prefix, float(time), agegrid_filename_suffix, agegrid_filename_ext),
            '-y {}'.format(time),
            '-n',
            # Use all feature types in proximity file (according to Dietmar)...
            #'-b',
            #'PassiveContinentalBoundary',
            '-x',
            '{}'.format(max_time),
            '-t',
            '1',
            '-i',
            '{}'.format(grid_spacing),
            #'-q',
            #str(proximity_threshold_kms),
            #'-d', # output distance with time
            '-j',
            '-w',
            '-c',
            str(1),
            '{}/distance_{}_{}'.format(output_dir, grid_spacing, time)])
    
    print('Time:', time)
    
    #print(' '.join(command_line))
    call_system_command(command_line)
    

    # Clamp the mean distance grids (and remove xy files).
    # Also rename the mean distance grids so that 'time' is at the end of the base filename -
    # this way we can import them as time-dependent raster into GPlates version 2.0 and earlier.
    #
    
    src_mean_distance_basename = '{}/distance_{}_{}_mean_distance'.format(output_dir, grid_spacing, time)
    dst_mean_distance_basename = '{}/mean_distance_{}d_{}'.format(output_dir, grid_spacing, time)
    
    src_mean_distance_xy = src_mean_distance_basename + '.xy'
    if os.access(src_mean_distance_xy, os.R_OK):
        os.remove(src_mean_distance_xy)
    
    src_mean_distance_grid = src_mean_distance_basename + '.nc'
    dst_mean_distance_grid = dst_mean_distance_basename + '.nc'
    
    if os.access(dst_mean_distance_grid, os.R_OK):
        os.remove(dst_mean_distance_grid)
    
    # Clamp mean distances.
    call_system_command(["gmt", "grdmath", "-fg", str(proximity_threshold_kms), src_mean_distance_grid, "MIN", "=", dst_mean_distance_grid])
    os.remove(src_mean_distance_grid)


# Wraps around 'generate_distance_grid()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def generate_distance_grid_parallel_pool_function(args):
    try:
        return generate_distance_grid(*args)
    except KeyboardInterrupt:
        pass


def low_priority():
    """ Set the priority of the process to below-normal."""

    import sys
    try:
        sys.getwindowsversion()
    except AttributeError:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        try:
            import psutil
        except ImportError:
            pass
        else:
            p = psutil.Process()
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    else:
        import os

        os.nice(1)


if __name__ == '__main__':
    
    print('Generating distance grids...')

    times = range(min_time, max_time + 1, time_step)
    #times = range(max_time, min_time - 1, -time_step) # Go backwards (can see results sooner).

    if use_all_cpus:
    
        # If 'use_all_cpus' is a bool (and therefore must be True) then use all available CPUs...
        if isinstance(use_all_cpus, bool):
            try:
                num_cpus = multiprocessing.cpu_count()
            except NotImplementedError:
                num_cpus = 1
        # else 'use_all_cpus' is a positive integer specifying the number of CPUs to use...
        elif isinstance(use_all_cpus, int) and use_all_cpus > 0:
            num_cpus = use_all_cpus
        else:
            raise TypeError('use_all_cpus: {} is neither a bool nor a positive integer'.format(use_all_cpus))
        
        # Split the workload across the CPUs.
        pool = multiprocessing.Pool(num_cpus, initializer=low_priority)
        pool_map_async_result = pool.map_async(
                generate_distance_grid_parallel_pool_function,
                (
                    (
                        time,
                    ) for time in range(min_time, max_time + 1, time_step)
                    #) for time in range(max_time, min_time - 1, -time_step) # Go backwards (can see results sooner).
                ),
                1) # chunksize

        # Apparently if we use pool.map_async instead of pool.map and then get the results
        # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
        # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        pool_map_async_result.get(99999)

    else:
        for time in times:
            generate_distance_grid(time)
