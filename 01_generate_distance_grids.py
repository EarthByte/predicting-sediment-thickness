
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
output_dir_base = '/Users/nickywright/PostDoc/Projects/STELLAR/paleobathymetry/paleobathymetry_traditional/TRUNK_2022_v2/sediment_thickness_D17'

if not os.path.exists(output_dir_base):
    os.makedirs(output_dir_base)


proximity_features_files = [
	'input_data/Global_EarthByte_GeeK07_COBLineSegments_2016_v4.gpmlz', # this is included in this repository
]

# --- lcoation of files on your computer
# DON'T FORGET TO UPDATE ocean_basin_proximity.py!
# --- agegrids
age_grid_dir = '/Users/nickywright/Data/Age/Muller2019-Young2019-Cao2020_Agegrids/Muller2019-Young2019-Cao2020_netCDF'   # change folder name if needed
age_grid_filename = 'Muller2019-Young2019-Cao2020_AgeGrid-'    # everything before 'time'
age_grid_filename_ext = 'nc'   # generally 'nc', but sometimes is 'grd'. Do not include the period

# --- topologies and other files
topology_dir = '/Users/nickywright/repos/usyd/EarthBytePlateMotionModel-ARCHIVE/Global_Model_WD_Internal_Release_2022_v2'
rotation_filenames = [
    '%s/Global_250-0Ma_Rotations.rot' % topology_dir,
    '%s/Global_410-250Ma_Rotations.rot' % topology_dir,
    '%s/Alps_Mesh_Rotations.rot' % topology_dir,
    '%s/Andes_Flat_Slabs_Rotations.rot' % topology_dir,
    '%s/Andes_Rotations.rot' % topology_dir,
    '%s/Australia_Antarctica_Mesh_Rotations.rot' % topology_dir,
    '%s/Australia_North_Zealandia_Rotations.rot' % topology_dir,
    '%s/Eurasia_Arabia_Mesh_Rotations.rot' % topology_dir,
    '%s/North_America_Flat_Slabs_Rotations.rot' % topology_dir,
    '%s/North_America_Mesh_Rotations.rot' % topology_dir,
    '%s/North_China_Mesh_Rotations.rot' % topology_dir,
    '%s/South_Atlantic_Rotations.rot' % topology_dir,
    '%s/South_China_DeformingModel.rot' % topology_dir,
    '%s/Southeast_Asia_Rotations.rot' % topology_dir,
]

topology_filenames = [
    '%s/Alps_Deforming_Mesh.gpml' % topology_dir,
    '%s/Alps_Mesh_Topologies.gpml' % topology_dir,
    '%s/America_Anyui_Deforming_Mesh.gpml' % topology_dir,
    '%s/America_Anyui_Mesh_Topologies.gpml' % topology_dir,
    '%s/Andes_Deforming_Mesh.gpml' % topology_dir,
    '%s/Andes_Flat_Slabs_Topologies.gpml' % topology_dir,
    '%s/Andes_Mesh_Topologies.gpml' % topology_dir,
    '%s/Arctic_Eurasia_Deforming_Mesh.gpml' % topology_dir,
    '%s/Arctic_Eurasia_Mesh_Topologies.gpml' % topology_dir,
    '%s/Australia_Antarctica_Deforming_Mesh.gpml' % topology_dir,
    '%s/Australia_Antarctica_Mesh_Topologies.gpml' % topology_dir,
    '%s/Australia_North_Zealandia_Deforming_Mesh.gpml' % topology_dir,
    '%s/Australia_North_Zealandia_Mesh_Topologies.gpml' % topology_dir,
    '%s/Baja_Deforming_Mesh.gpml' % topology_dir,
    '%s/Coral_Sea_Deforming_Mesh.gpml' % topology_dir,
    '%s/Coral_Sea_Topologies.gpml' % topology_dir,
    '%s/East_African_Rift_Deforming_Mesh_and_Topologies.gpml' % topology_dir,
    '%s/East-West_Gondwana_Deforming_Mesh_and_Topologies.gpml' % topology_dir,
    '%s/Ellesmere_Deforming_Mesh.gpml' % topology_dir,
    '%s/Eurasia_Arabia_Deforming_Mesh_and_Topologies.gpml' % topology_dir,
    '%s/Global_Mesozoic-Cenozoic_PlateBoundaries.gpml' % topology_dir,
    '%s/Global_Paleozoic_PlateBoundaries.gpml' % topology_dir,
    '%s/Greater_India_Deforming_Mesh.gpml' % topology_dir,
    '%s/Greater_India_Mesh_Topologies.gpml' % topology_dir,
    '%s/Inactive_Meshes_and_Topologies.gpml' % topology_dir,
    '%s/North_America_Mesh_Topologies.gpml' % topology_dir,
    '%s/North_Atlantic_Deforming_Mesh.gpml' % topology_dir,
    '%s/North_Atlantic_Mesh_Topologies.gpml' % topology_dir,
    '%s/North_China_Mesh_Topologies.gpml' % topology_dir,
    '%s/Northern_Andes_Deforming_Mesh.gpml' % topology_dir,
    '%s/Northern_Andes_Mesh_Topologies.gpml' % topology_dir,
    '%s/Papua_New_Guinea_Deforming_Meshes.gpml' % topology_dir,
    '%s/Papua_New_Guinea_Mesh_Topologies.gpml' % topology_dir,
    '%s/Scotia_Deforming_Mesh_and_Topologies.gpml' % topology_dir,
    '%s/Siberia_Eurasia_Deforming_Mesh.gpml' % topology_dir,
    '%s/Siberia_Eurasia_Mesh_Topologies.gpml' % topology_dir,
    '%s/South_Atlantic_Deforming_Mesh.gpml' % topology_dir,
    '%s/South_Atlantic_Mesh_Topologies.gpml' % topology_dir,
    '%s/South_China_Mesh_Topologies.gpml' % topology_dir,
    '%s/South_China_DeformingElements.gpml' % topology_dir,
    '%s/South_Zealandia_Deforming_Mesh.gpml' % topology_dir,
    '%s/South_Zealandia_Mesh_Topologies.gpml' % topology_dir,
    '%s/Southeast_Asia_Deforming_Mesh.gpml' % topology_dir,
    '%s/Southeast_Asia_Mesh_Topologies.gpml' % topology_dir,
    '%s/West_Antarctic_Zealandia_Deforming_Mesh.gpml' % topology_dir,
    '%s/West_Antarctica_Zealandia_Mesh_Topologies.gpml' % topology_dir,
    '%s/Western_North_America_Deforming_Mesh.gpml' % topology_dir,
    '%s/Western_Tethys_Deforming_Mesh.gpml' % topology_dir,
    '%s/Western_Tethys_Tectonic_Boundary_Topologies.gpml' % topology_dir]


# --- set times and spacing
grid_spacing = 0.2

min_time = 0
max_time = 250
time_step = 1

proximity_threshold_kms = 3000

output_dir = '%s/distances_%sd' % (output_dir_base, grid_spacing)

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
        num_cpus = multiprocessing.cpu_count() - 2
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
