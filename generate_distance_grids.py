from __future__ import print_function
from call_system_command import call_system_command
import math
import multiprocessing
import os
import sys

output_dir = 'distances_1d'

topology_dir = 'E:/Users/John/Downloads/GPlates/data/vector/Muller_etal_AREPS_Supplement'
#topology_dir = 'D:/Users/john/Downloads/gplates/data/PlateModels/Muller_etal_AREPS_Supplement'
age_grid_dir = 'E:/Users/John/Downloads/GPlates/data/rasters/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'
#age_grid_dir = 'D:/Users/john/Downloads/gplates/data/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'

proximity_features_file = 'input_data/Global_EarthByte_GeeK07_COBLineSegments_2016_v4.gpmlz'

grid_spacing = 1.0

min_time = 0
max_time = 230
time_step = 1

proximity_threshold_kms = 3000


def generate_distance_grid(time):
    
    command_line = [
            'python',
            'ocean_basin_proximity.py',
            '-r',
            '{0}/Global_EarthByte_230-0Ma_GK07_AREPS.rot'.format(topology_dir),
            '-m',
            proximity_features_file,
            '-s',
            '{0}/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml'.format(topology_dir),
            '{0}/Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml'.format(topology_dir),
            '-g',
            '{0}/EarthByte_AREPS_Muller_etal_2016_AgeGrid-{1}.nc'.format(age_grid_dir, time),
            '-y {0}'.format(time),
            '-n',
            # Use all feature types in proximity file (according to Dietmar)...
            #'-b',
            #'PassiveContinentalBoundary',
            '-x',
            '230',
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
            '{0}/distance_{1}_{2}'.format(output_dir, grid_spacing, time)]
    
    print('Time:', time)
    
    call_system_command(command_line)
    
    #
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
