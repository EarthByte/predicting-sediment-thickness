
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
    
    # Sample the sediment thickness raster at the subduction points.
    sediment_thicknesses = raster_query.query_raster_at_points(
            sediment_thickness_grid_filename,
            subduction_points,
            search_radius_radians,
            smoothing_radius_radians)
    
    subducting_lon_lat_thickness_velocity_volume_list = []
    
    # Iterate over subduction points/thicknesses and calculate statistics (including subducting volume).
    weighted_mean_subducting_sed_thickness = 0.0
    weighted_second_moment_subducting_sed_thickness = 0.0
    total_subducting_length_metres = 0.0
    total_subducting_sediment_volume_metres_3_per_year = 0.0
    for subduction_point_index, sediment_thickness in enumerate(sediment_thicknesses):
        if math.isnan(sediment_thickness):
            continue
        
        subduction_convergence_item = subduction_convergence_data[subduction_point_index]
        
        subducting_lon = subduction_convergence_item[0]
        subducting_lat = subduction_convergence_item[1]
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
        
        weighted_mean_subducting_sed_thickness += subducting_length_metres * sediment_thickness
        weighted_second_moment_subducting_sed_thickness += subducting_length_metres * sediment_thickness * sediment_thickness
        
        total_subducting_length_metres += subducting_length_metres
        convergence_normal_velocity_metres_per_year = (
            # 1e-2 converts cm/y to m/y...
            1e-2 * math.fabs(convergence_velocity_magnitude_cm_per_yr) *
            # Negative convergence handled by cos(obliquity_angle)...
            math.cos(math.radians(convergence_obliquity_degrees)))
        subducting_sediment_volume_metres_3_per_year = (
                sediment_thickness * subducting_length_metres * convergence_normal_velocity_metres_per_year)
        total_subducting_sediment_volume_metres_3_per_year += subducting_sediment_volume_metres_3_per_year
        
        subducting_lon_lat_thickness_velocity_volume_list.append((
                subducting_lon,
                subducting_lat,
                sediment_thickness,
                subducting_sediment_volume_metres_3_per_year / subducting_length_metres,
                # cms/year ...
                1e2 * convergence_normal_velocity_metres_per_year))
    
    # mean = M = sum(Ci * Xi) / sum(Ci)
    # std_dev  = sqrt[sum(Ci * (Xi - M)^2) / sum(Ci)]
    #          = sqrt[(sum(Ci * Xi^2) - 2 * M * sum(Ci * Xi) + M^2 * sum(Ci)) / sum(Ci)]
    #          = sqrt[(sum(Ci * Xi^2) - 2 * M * M * sum(Ci) + M^2 * sum(Ci)) / sum(Ci)]
    #          = sqrt[(sum(Ci * Xi^2) - M^2 * sum(Ci)) / sum(Ci)]
    #          = sqrt[(sum(Ci * Xi^2) / sum(Ci) - M^2]
    mean_sed_thickness = weighted_mean_subducting_sed_thickness / total_subducting_length_metres
    variance_sed_thickness = (
            (weighted_second_moment_subducting_sed_thickness / total_subducting_length_metres) -
            mean_sed_thickness * mean_sed_thickness)
    std_dev_sed_thickness = math.sqrt(variance_sed_thickness) if variance_sed_thickness > 0.0 else 0.0
    
    return (time,
            subducting_lon_lat_thickness_velocity_volume_list,
            mean_sed_thickness,
            std_dev_sed_thickness,
            total_subducting_length_metres,
            total_subducting_sediment_volume_metres_3_per_year)


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
    
    # Print the header of the statistics table.
    print('{0:<10} {1:<40} {2:<40} {3:<40} {4:<40}'.format(
            'Age(Ma)',
            'Mean subducting thickness (m)',
            'Std dev subducting thickness (m)',
            'Total subducting volume (m^3/y)',
            'Subducting volume per unit metre (m^2/y)'))
    
    subducting_thickness_features = []
    
    # Iterate over the statistics for each time (each time comes from a separate multiprocessing pool).
    for (time,
        subducting_lon_lat_thickness_velocity_volume_list,
        mean_sed_thickness,
        std_dev_sed_thickness,
        total_subducting_length_metres,
        total_subducting_sediment_volume_metres_3_per_year) in sorted(subduction_datas):
        
        # Print the statistics for the current time.
        print('{0:<10.0f} {1:<40.2f} {2:<40.2f} {3:<40.2f} {4:<40.2f}'.format(
                time,
                mean_sed_thickness,
                std_dev_sed_thickness,
                total_subducting_sediment_volume_metres_3_per_year,
                total_subducting_sediment_volume_metres_3_per_year / total_subducting_length_metres))
        
        # Gather the subducting thickness points.
        subducting_points = []
        subducting_sed_thicknesses = []
        subducting_sediment_volumes_metres_3_per_year_per_metre = []
        convergence_normal_velocities_cms_per_year = []
        for (
            lon,
            lat,
            sed_thickness,
            subducting_sediment_volume_metres_3_per_year_per_metre,
            convergence_normal_velocity_cms_per_year) in subducting_lon_lat_thickness_velocity_volume_list:
            
            subducting_points.append(pygplates.PointOnSphere(lat, lon))
            subducting_sed_thicknesses.append(sed_thickness)
            subducting_sediment_volumes_metres_3_per_year_per_metre.append(subducting_sediment_volume_metres_3_per_year_per_metre)
            convergence_normal_velocities_cms_per_year.append(convergence_normal_velocity_cms_per_year)
        
        # Create a scalar coverage feature to display sediment thicknesses in GPlates.
        subducting_thickness_feature = pygplates.Feature()
        subducting_thickness_feature.set_geometry((
                pygplates.MultiPointOnSphere(subducting_points),
                {pygplates.ScalarType.create_gpml('subducting_sed_thick') : subducting_sed_thicknesses,
                pygplates.ScalarType.create_gpml('sed_volume_m_3_per_year_per_m') : subducting_sediment_volumes_metres_3_per_year_per_metre,
                pygplates.ScalarType.create_gpml('conv_normal_vel_cms_year') : convergence_normal_velocities_cms_per_year}))
        # Only want to display this feature at 'time' Ma.
        subducting_thickness_feature.set_valid_time(time + 0.5, time - 0.5)
        
        subducting_thickness_features.append(subducting_thickness_feature)
    
    pygplates.FeatureCollection(subducting_thickness_features).write('subducting_thicknesses.gpmlz')
    
    sys.exit(0)