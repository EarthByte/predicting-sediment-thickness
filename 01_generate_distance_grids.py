
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

# For each distance grid do not reconstruct ocean points earlier than 'max_topological_reconstruction_time'
# (each ocean point is reconstructed back to its age grid value or this value, whichever is smaller).
# This limit can be set to the earliest (max) reconstruction time of the topological model.
# If it's 'None' then only the age grid limits how far back each point is reconstructed.
max_topological_reconstruction_time = 1000  # can be None to just use age grid as the limit

# Distances in the mean distance grids are clamped to this value (in kms).
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
def generate_distance_grid(times):
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

    # Age grid filenames and paleo times.
    age_grid_filenames_and_times = []
    for time in times:
        age_grid_filename = '{}/{}{:.1f}{}.{}'.format(agegrid_dir, agegrid_filename_prefix, time, agegrid_filename_suffix, agegrid_filename_ext)
        age_grid_filenames_and_times.append(age_grid_filename)
        age_grid_filenames_and_times.append(str(time))
    command_line.append('-g')
    command_line.extend(age_grid_filenames_and_times)

    # Use all feature types in proximity file (according to Dietmar)...
    #command_line.extend(['-b', 'PassiveContinentalBoundary'])

    # Time increment is 1 Myr (this is for topological reconstruction of the ocean points).
    command_line.extend(['-t', '1'])
    
    # If limiting the max topological reconstruction time.
    if max_topological_reconstruction_time is not None:
        command_line.extend(['-x', '{}'.format(max_topological_reconstruction_time)])
    
    # Grid spacing.
    command_line.extend(['-i', '{}'.format(grid_spacing)])

    # Proximity threshold - but we're not using here - instead clamping the mean distance grids explicitly below.
    #command_line.extend(['-q', str(proximity_threshold_kms)])

    # Don't output distance grids for all reconstruction times.
    # Only outputting a single "mean" (over all reconstruction times) distance grid.
    #command_line.append('-d')

    # Output a "mean" (over all reconstruction times) distance grid.
    command_line.append('-j')

    # Generate a grd (".nc") file for each xyz file.
    command_line.append('-w')
    
    # Number of cores.
    command_line.extend(['-c', '1'])

    # Distance grids output filename prefix.
    command_line.append('{}/distance_{}'.format(output_dir, grid_spacing))
    
    print('Times:', list(times))
    
    #print(' '.join(command_line))
    call_system_command(command_line)
    

    # Clamp the mean distance grids (and remove xy files).
    # Also rename the mean distance grids so that 'time' is at the end of the base filename -
    # this way we can import them as time-dependent raster into GPlates version 2.0 and earlier.
    #
    
    for time in  times:
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

    #import time as time_prof

    print('Generating distance grids...')
    
    #tprof_start = time_prof.perf_counter()

    times = range(min_time, max_time + 1, time_step)
    #times = range(max_time, min_time - 1, -time_step) # Go backwards (can see results sooner).

    # Give each task a reasonable number of age grids (times) to process - if there's not enough times per task then we'll
    # spend too much time resolving/reconstructing proximity features (and generating shortest path obstacle grids) -
    # which needs to be repeated for each task (group of times) - this is because each age grid involves
    # reconstructing all its ocean points back in time until they disappear (at mid-ocean ridge) and so there's
    # a lot of overlap in time across the age grids (where calculations can be shared within a single process).
    min_num_age_grids_per_call = 10

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
        
        # Divide the input times (age grid paleo times) into sub-lists (each to be run on a separate process).
        #
        # We could reduce the number of tasks to the number of CPUs (ie, increase number of times per task).
        # However some tasks might finish sooner than others leaving some CPUs under utilised. Conversely, we don't
        # want the number of tasks to be too high otherwise the calculation sharing (mentioned above) is reduced.
        # So we double the number of tasks (twice number of CPUs) as a good compromise.
        num_tasks = 2 * num_cpus
        
        # Create the pool sub-lists of times.
        times_sub_lists = []
        num_times = len(times)
        num_times_per_task = (num_times + num_tasks - 1) // num_tasks
        if num_times_per_task < min_num_age_grids_per_call:
            num_times_per_task = min_num_age_grids_per_call
        task_start_time_index = 0
        while task_start_time_index < num_times:
            # Pass consecutive times into each process since the distance calculations are more efficient that way.
            times_sub_lists.append(times[task_start_time_index : task_start_time_index + num_times_per_task])
            task_start_time_index += num_times_per_task
        
        try:
            # Split the workload across the CPUs.
            pool = multiprocessing.Pool(num_cpus, initializer=low_priority)
            pool_map_async_result = pool.map_async(
                    generate_distance_grid_parallel_pool_function,
                    (
                        (
                            times_sub_list,
                        ) for times_sub_list in times_sub_lists
                    ),
                    1) # chunksize
            
            # Apparently if we use pool.map_async instead of pool.map and then get the results
            # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
            # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
            pool_map_async_result.get(999999)
        except KeyboardInterrupt:
            # Note: 'finally' block below gets executed before returning.
            pass
        finally:
            pool.close()
            pool.join()

    else:
        num_times = len(times)
        start_time_index = 0
        while start_time_index < num_times:
            # Pass 'min_num_age_grids_per_call' consecutive times into each process since the distance calculations are more efficient that way.
            generate_distance_grid(times[start_time_index : start_time_index + min_num_age_grids_per_call])
            start_time_index += min_num_age_grids_per_call
    
    #tprof_end = time_prof.perf_counter()
    #print(f"Total time: {tprof_end - tprof_start:.2f} seconds")
