from __future__ import print_function
import os
import sys
from call_system_command import call_system_command

output_dir = 'distances_1d'

topology_dir = 'E:/Users/John/Downloads/GPlates/data/vector/Muller_etal_AREPS_Supplement'
#topology_dir = 'D:/Users/john/Downloads/gplates/data/PlateModels/Muller_etal_AREPS_Supplement'
age_grid_dir = 'E:/Users/John/Downloads/GPlates/data/rasters/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'
#age_grid_dir = 'D:/Users/john/Downloads/gplates/data/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'

proximity_features_file = 'input_data/Global_EarthByte_GeeK07_COBLineSegments_2016_v4.gpmlz'

grid_spacing = 1

for time in range(0, 231):
    
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
            '{0}/distance_{1}_{2}'.format(output_dir, grid_spacing, time)]
    
    print('Time:', time)
    
    # Capture stderr since pygplates generates a bunch of topology warning messages which we just ignore.
    stderr = call_system_command(command_line, return_stderr=True)
    
    # Rename the mean distance files so that 'time' is at the end of the base filename -
    # this way we can import them as time-dependent raster into GPlates.
    for ext in ('xy', 'grd'):
        
        src_mean_distance = '{0}/distance_{1}_{2}_mean_distance.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_mean_distance = '{0}/mean_distance_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        #print(src_mean_distance, dst_mean_distance)
        
        if os.access(dst_mean_distance, os.R_OK):
            os.remove(dst_mean_distance)
        os.rename(src_mean_distance, dst_mean_distance)
