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
            '{0}/agegrid_{1}.nc'.format(age_grid_dir, time),
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
            '-j',
            '-w',
            '-c',
            str(1),
            '{0}/distance_{1}_{2}'.format(output_dir, grid_spacing, time)]
    
    print('Time:', time)
    
    call_system_command(command_line)
    
    # Rename the mean distance files so that 'time' is at the end of the base filename -
    # this way we can import them as time-dependent raster into GPlates.
    for ext in ('xy', 'grd'):
        
        src_mean_distance = '{0}/distance_{1}_{2}_mean_distance.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_mean_distance = '{0}/mean_distance_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        #print(src_mean_distance, dst_mean_distance)
        
        if os.access(dst_mean_distance, os.R_OK):
            os.remove(dst_mean_distance)
        os.rename(src_mean_distance, dst_mean_distance)


# Wraps around 'generate_distance_grid()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def generate_distance_grid_parallel_pool_function(args):
    try:
        return generate_distance_grid(*args)
    except KeyboardInterrupt:
        pass



if __name__ == '__main__':
    
    #generate_distance_grid(197)
    #sys.exit(0)
    
    try:
        num_cpus = multiprocessing.cpu_count()
    except NotImplementedError:
        num_cpus = 1
    
    #num_cpus = 11
    
    print('Generating distance grids...')
    
    # Split the workload across the CPUs.
    pool = multiprocessing.Pool(num_cpus)
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
