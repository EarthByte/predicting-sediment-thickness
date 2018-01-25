from __future__ import print_function
import math
import multiprocessing
import os
import pygplates
import raster_query
import subduction_convergence
import sys

rotation_filename = 'E:/Users/John/Downloads/GPlates/data/vector/Muller_etal_AREPS_Supplement/Global_EarthByte_230-0Ma_GK07_AREPS.rot'

topology_dir = 'E:/Users/John/Downloads/GPlates/data/vector/Muller_etal_AREPS_Supplement'
#topology_dir = 'D:/Users/john/Downloads/gplates/data/PlateModels/Muller_etal_AREPS_Supplement'
sediment_thickness_grid_dir = 'E:/Users/John/Downloads/GPlates/data/PythonWorkflows/SedimentationRate/sedimentation_output/predicted_thickness'

tessellation_threshold_radians = math.radians(0.5)
search_radius_radians = math.radians(20.0)
smoothing_radius_radians = None

min_time = 0
max_time = 230
time_step = 1


def calc_subducting_sediment_volume(time):

    rotation_model = pygplates.RotationModel(rotation_filename)
    topology_features = pygplates.FeaturesFunctionArgument(
            (topology_dir + '/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml',
            topology_dir + '/Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml')).get_features()

    sediment_thickness_grid_filename = sediment_thickness_grid_dir + '/sed_thick_0.2d_{0}.grd'.format(time)
    
    subduction_convergence_data = subduction_convergence.subduction_convergence(
                rotation_model,
                topology_features,
                tessellation_threshold_radians,
                time)
    
    subduction_points = [pygplates.PointOnSphere(data[1], data[0]) for data in subduction_convergence_data]
    
    sediment_thicknesses = raster_query.query_raster_at_points(
            sediment_thickness_grid_filename,
            subduction_points,
            search_radius_radians,
            smoothing_radius_radians)
    
    #sediment_thickness_data = raster_query.query_raster_with_resolved_topologies(
    #        sediment_thickness_grid_filename,
    #        rotation_model,
    #        time,
    #        topology_features,
    #        tessellation_threshold_radians,
    #        search_radius_radians,
    #        smoothing_radius_radians,
    #        [pygplates.FeatureType.gpml_subduction_zone])
    
    total_subducting_length_metres = 0.0
    total_subducting_sediment_volume_metres_3_per_year = 0.0
    #count = 0
    for subduction_point_index, sediment_thickness in enumerate(sediment_thicknesses):
        if math.isnan(sediment_thickness):
            #count += 1
            continue
        
        subduction_convergence_item = subduction_convergence_data[subduction_point_index]
        
        #lon = subduction_convergence_item[0]
        #lat = subduction_convergence_item[1]
        convergence_velocity_magnitude_cm_per_yr = subduction_convergence_item[2]
        convergence_obliquity_degrees = subduction_convergence_item[3]
        #absolute_velocity_magnitude = subduction_convergence_item[4]
        #absolute_obliquity_degrees = subduction_convergence_item[5]
        subducting_length_degrees = subduction_convergence_item[6]
        #subducting_arc_normal_azimuth = subduction_convergence_item[7]
        #subducting_plate_id = subduction_convergence_item[8]
        #overriding_plate_id = subduction_convergence_item[9]
        
        subducting_length_metres = (
            math.radians(subducting_length_degrees) * 1e3 * pygplates.Earth.mean_radius_in_kms)
        
        total_subducting_length_metres += subducting_length_metres
        
        convergence_normal_velocity_metres_per_year = (
            # 1e-2 converts cm/y to m/y...
            1e-2 * math.fabs(convergence_velocity_magnitude_cm_per_yr) *
            # Negative convergence handled by cos(obliquity_angle)...
            math.cos(math.radians(convergence_obliquity_degrees)))
        
        total_subducting_sediment_volume_metres_3_per_year += (
                sediment_thickness * subducting_length_metres * convergence_normal_velocity_metres_per_year)
    
    #print('Count', count, len(sediment_thicknesses))
    
    #mean_sediment_thickness = sum(sediment_thicknesses) / len(sediment_thicknesses)
    
    return time, total_subducting_length_metres, total_subducting_sediment_volume_metres_3_per_year


# Wraps around 'calc_subducting_sediment_volume()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def calc_subducting_sediment_volume_parallel_pool_function(args):
    try:
        return calc_subducting_sediment_volume(*args)
    except KeyboardInterrupt:
        pass



if __name__ == '__main__':
    
    #calc_subducting_sediment_volume(0)
    #sys.exit(0)
    
    try:
        num_cpus = multiprocessing.cpu_count()
        #num_cpus = 10
    except NotImplementedError:
        num_cpus = 1
    
    # Split the workload across the CPUs.
    pool = multiprocessing.Pool(num_cpus)
    pool_map_async_result = pool.map_async(
            calc_subducting_sediment_volume_parallel_pool_function,
            (
                (
                    time,
                ) for time in range(min_time, max_time + 1, time_step)
            ),
            1) # chunksize

    # Apparently if we use pool.map_async instead of pool.map and then get the results
    # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
    # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    try:
        subduction_datas = pool_map_async_result.get(99999)
    except KeyboardInterrupt:
        sys.exit(1)
    
    print('{0:<10} {1:<40} {2:<40}'.format('Age(Ma)', 'Total subducting volume (m^3/y)', 'Subducting volume per unit metre (m^2/y)'))
    
    for time, total_subducting_length_metres, total_subducting_sediment_volume_metres_3_per_year in sorted(subduction_datas):
        print('{0:<10.0f} {1:<40.2f} {2:<40.2f}'.format(
                time,
                total_subducting_sediment_volume_metres_3_per_year,
                total_subducting_sediment_volume_metres_3_per_year / total_subducting_length_metres))
    
    sys.exit(0)