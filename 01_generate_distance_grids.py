
from call_system_command import call_system_command
import math
import multiprocessing
import os, shutil
import sys

""" This script creates grids (netcdfs) of the mean distance to passive margins through time


Requirements & Inputs:
    - Python  
    - Python scripts: ocean_basin_proximity.py, 
    - GMT 5 (or later)
    - Files associated with a tectonic model: agegrids, rotation file, plate boundaries, topologies

Outputs:
    - directory named 'distances_1d', with mean distance grids (in metres) through time.

2020-02-14: Added comments, removed hardcoded names (for rotation file, etc) from definitions
"""

# ----- set directories and filenames

output_dir = 'distances_1d'


proximity_features_files = [
	'input_data/Global_EarthByte_GeeK07_COBLineSegments_2016_v4.gpmlz', # this is included in this repository
]

# --- lcoation of files on your computer
data_dir = '/home/michael/workspace/predicting-sediment-thickness/Muller_etal_2016_AREPS'

# --- agegrids
age_grid_dir = '%s/Muller_etal_2016_AREPS_Agegrids/Muller_etal_2016_AREPS_Agegrids_v1.17/Muller_etal_2016_AREPS_v1.17_netCDF' % data_dir   # change folder name if needed
age_grid_filename = 'Muller_etal_2016_AREPS_v1.17_AgeGrid-'    # everything before 'time'
age_grid_filename_ext = 'nc'   # generally 'nc', but sometimes is 'grd'. Do not include the period

# --- topologies and other files
topology_dir = '%s' % data_dir
rotation_filenames = [
	'%s/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS.rot' % topology_dir,
]
topology_filenames = [
	'%s/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml' % topology_dir,
	'%s/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml' % topology_dir,
]


# --- set times and spacing
grid_spacing = 1.0

min_time = 0
max_time = 230
time_step = 1

proximity_threshold_kms = 3000

# -----
if not os.path.exists(output_dir):
    print('%s does not exist, creating now... ' % output_dir)
    os.mkdir(output_dir)

# ----- 
def generate_distance_grid(time):
    py_cmd='python3'
    if shutil.which('python3') is None:
        py_cmd = 'python'
    
    command_line = [py_cmd, 'ocean_basin_proximity.py']
    command_line.extend(['-r'])
    command_line.extend('{0}'.format(rotation_filename) for rotation_filename in rotation_filenames)
    command_line.extend(['-m'])
    command_line.extend('{0}'.format(proximity_features_file) for proximity_features_file in proximity_features_files)
    command_line.extend(['-s'])
    command_line.extend('{0}'.format(topology_filename) for topology_filename in topology_filenames)
    command_line.extend([
            '-g',
            '{0}/{1}{2}.{3}'.format(age_grid_dir, age_grid_filename, time, age_grid_filename_ext),
            '-y {0}'.format(time),
            '-n',
            # Use all feature types in proximity file (according to Dietmar)...
            #'-b',
            #'PassiveContinentalBoundary',
            '-x',
            '{0}'.format(max_time),
            '-t',
            '1',
            '-i',
            '{0}'.format(grid_spacing),
            #'-q',
            #str(proximity_threshold_kms),
            #'-d', # output distance with time
            '-j',
            '-w',
            '-c',
            str(1),
            '{0}/distance_{1}_{2}'.format(output_dir, grid_spacing, time)])
    
    print('Time:', time)
    
    #print(' '.join(command_line))
    call_system_command(command_line)
    

    # Clamp the mean distance grids (and remove xy files).
    # Also rename the mean distance grids so that 'time' is at the end of the base filename -
    # this way we can import them as time-dependent raster into GPlates version 2.0 and earlier.
    #
    
    src_mean_distance_basename = '{0}/distance_{1}_{2}_mean_distance'.format(output_dir, grid_spacing, time)
    dst_mean_distance_basename = '{0}/mean_distance_{1}d_{2}'.format(output_dir, grid_spacing, time)
    
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
        import psutil
        
        p = psutil.Process()
        p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    else:
        import os

        os.nice(1)


if __name__ == '__main__':
    
    try:
        num_cpus = multiprocessing.cpu_count()
    except NotImplementedError:
        num_cpus = 1
    
    print('Generating distance grids...')
    
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
