from __future__ import print_function
from call_system_command import call_system_command
import multiprocessing
import os
import sys

age_grid_dir = 'E:/Users/John/Downloads/GPlates/data/rasters/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'
#age_grid_dir = 'D:/Users/john/Downloads/gplates/data/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma'
distance_grid_dir = 'distances_1d'

grid_spacing = 0.2

min_time = 0
max_time = 230
time_step = 1


def generate_predicted_sedimentation_grid(
        time,
        predict_sedimentation_script,
        scale_sedimentation_rate,
        mean_age,
        mean_distance,
        variance_age,
        variance_distance,
        max_age,
        max_distance,
        age_distance_polynomial_coefficients,
        output_dir):
    
    command_line = [
            'python',
            predict_sedimentation_script,
            '-d',
            '{0}/mean_distance_1d_{1}.grd'.format(distance_grid_dir, time),
            '-g',
            '{0}/agegrid_{1}.nc'.format(age_grid_dir, time),
            '-i',
            str(grid_spacing),
            '-w',
            '-m',
            str(mean_age),
            str(mean_distance),
            '-v',
            str(variance_age),
            str(variance_distance),
            '-x',
            str(max_age),
            str(max_distance),
            '-f']
    command_line.extend(str(coeff) for coeff in age_distance_polynomial_coefficients)
    command_line.extend([
            '-p',
            '0.63', # Surface porosity (shale)
            '-c',
            '5.71e-4', # Porosity decay (shale)
            '-s',
            str(scale_sedimentation_rate),
            '--',
            '{0}/sed_{1}_{2}'.format(output_dir, grid_spacing, time)])
    
    #print('Time:', time)
    #print(command_line)
    
    # Execute the command.
    call_system_command(command_line)
    
    # Rename the average sedimentation rate and sediment thicknesses files so that 'time' is at the
    # end of the base filename - this way we can import them as time-dependent raster into GPlates.
    for ext in ('xy', 'grd'):
        
        src_sed_rate = '{0}/sed_{1}_{2}_sed_rate.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_sed_rate = '{0}/sed_rate_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        if os.access(dst_sed_rate, os.R_OK):
            os.remove(dst_sed_rate)
        os.rename(src_sed_rate, dst_sed_rate)
        
        src_sed_thick = '{0}/sed_{1}_{2}_sed_thick.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_sed_thick = '{0}/sed_thick_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        if os.access(dst_sed_thick, os.R_OK):
            os.remove(dst_sed_thick)
        os.rename(src_sed_thick, dst_sed_thick)


# Wraps around 'generate_predicted_sedimentation_grid()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def generate_predicted_sedimentation_grid_parallel_pool_function(args):
    try:
        return generate_predicted_sedimentation_grid(*args)
    except KeyboardInterrupt:
        pass



if __name__ == '__main__':
    
    try:
        num_cpus = multiprocessing.cpu_count()
        num_cpus = 10
    except NotImplementedError:
        num_cpus = 1

    # Machine learning training parameters.
    # These come from the "sediment_rate_decomp" IPython notebook (using sklearn Python module).

    #
    # Predict sedimentation *rate* results in "sediment_rate_decomp_v5.ipynb" from:
    #     
    #     _, regressor_trained_no_river =  geo_preprocess3.regression(data=data_train, 
    #                                                                     regressor=regressor_no_river, 
    #                                                                     n_splits=n_splits,
    #                                                                     lon_ind=lon_ind, 
    #                                                                     lat_ind=lat_ind, 
    #                                                                     y_ind=-1,
    #                                                                     logy=True)
    #     
    #     print('Mean of age and distance:', regressor_trained_no_river.named_steps['stand'].mean_)
    #     print('Variance of age and distance:', regressor_trained_no_river.named_steps['stand'].var_)
    #     print('Polynomial coefficients:', regressor_trained_no_river.named_steps['linear'].coef_)
    #     print('Polynomial intercept:', regressor_trained_no_river.named_steps['linear'].intercept_)
    #     print('Polynomial feature names:', regressor_trained_no_river.named_steps['poly'].get_feature_names())
    #
    #
    # predict_sedimentation_script = 'predict_sedimentation_rate.py'
    # scale_sedimentation_rate = 10.0  # Scale predicted rate (cm/Ky) to (m/My).
    # mean_age = 57.08516053
    # mean_distance = 2004.37249998
    # variance_age = 1.57169637e+03
    # variance_distance = 2.43160312e+06
    # age_distance_polynomial_coefficients = [
    #         0.0 , -0.52275531, -0.51023915,  0.34082993, -0.08491046, 0.5764176 , -0.0704285 , -0.01460767,  0.1403967 , -0.24019863]

    #
    # Predict sedimentation *rate* results in "sediment_rate_decomp_v5.ipynb" from:
    #
    #     dataq = data[:, [lon_ind, lat_ind, age_ind, passive_dis_ind, sedrate_ind]]
    #     geo_preprocess3.two_feature_analysis(dataq, regressor, 2, 3, 'age', 
    #                                          'distance to passive margin', 'predicted log sedrate',
    #                                          query_size=20)
    #     
    #     print('Mean of age and distance:', regressor.named_steps['stand'].mean_)
    #     print('Variance of age and distance:', regressor.named_steps['stand'].var_)
    #     print('Polynomial coefficients:', regressor.named_steps['linear'].coef_)
    #     print('Polynomial intercept:', regressor.named_steps['linear'].intercept_)
    #     print('Polynomial feature names:', regressor.named_steps['poly'].get_feature_names())
    #
    #
    predict_sedimentation_script = 'predict_sedimentation_rate.py'
    #scale_sedimentation_rate = 1.0  # Keep predicted rate in (cm/Ky).
    scale_sedimentation_rate = 10.0  # Scale predicted rate (cm/Ky) to (m/My).
    mean_age = 60.19000622
    mean_distance = 2058.25803042
    variance_age = 1.89421205e+03
    variance_distance = 2.96811661e+06
    max_age = 196.88598633
    max_distance = 6855.48388672
    age_distance_polynomial_coefficients = [
            -1.3607332291017404 , -0.37400809, -0.5160163 ,  0.42959   , -0.07323361, 0.60509602, -0.12621648, -0.01720338,  0.14899377, -0.29455261]
    
    output_dir = 'sedimentation_output/predicted_rate'
    
    print('Generating predicted sedimentation rate grids...')
    
    # Split the workload across the CPUs.
    pool = multiprocessing.Pool(num_cpus)
    pool_map_async_result = pool.map_async(
            generate_predicted_sedimentation_grid_parallel_pool_function,
            (
                (
                    time,
                    predict_sedimentation_script,
                    scale_sedimentation_rate,
                    mean_age,
                    mean_distance,
                    variance_age,
                    variance_distance,
                    max_age,
                    max_distance,
                    age_distance_polynomial_coefficients,
                    output_dir
                ) for time in range(min_time, max_time + 1, time_step)
            ),
            1) # chunksize

    # Apparently if we use pool.map_async instead of pool.map and then get the results
    # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
    # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    pool_map_async_result.get(99999)
    
    #
    # Predict sediment *thickness* results in "sediment_thick_v5.ipynb" from:
    #
    #     dataq = data[:, [lon_ind, lat_ind, age_ind, passive_dis_ind, sedthick_ind]]
    #     geo_preprocess3.two_feature_analysis(dataq, regressor, 2, 3, 'age', 
    #                                          'distance to passive margin', 'predicted log thickness',
    #                                          query_size=20)
    #     
    #     print('Mean of age and distance:', regressor.named_steps['stand'].mean_)
    #     print('Variance of age and distance:', regressor.named_steps['stand'].var_)
    #     print('Max of age and distance:', np.max(data[:, [age_ind, passive_dis_ind]], axis=0))
    #     print('Polynomial coefficients:', regressor.named_steps['linear'].coef_)
    #     print('Polynomial intercept:', regressor.named_steps['linear'].intercept_)
    #     print('Polynomial feature names:', regressor.named_steps['poly'].get_feature_names())
    #
    #
    predict_sedimentation_script = 'predict_sediment_thickness.py'
    scale_sedimentation_rate = 1.0  # No scaling - we're calculating rate (m/My) from thickness (m) and age (My).
    mean_age = 60.19000622
    mean_distance = 2058.25803042
    variance_age = 1.89421205e+03
    variance_distance = 2.96811661e+06
    max_age = 196.88598633
    max_distance = 6855.48388672
    age_distance_polynomial_coefficients = [
            5.0517840448654425, 0.35892409, -0.39697653, -0.22518593, -0.03911694, 0.57485024,  0.06178908, -0.04514372,  0.17936484, -0.29459847]
    
    output_dir = 'sedimentation_output/predicted_thickness'

    print('Generating predicted sediment thickness grids...')
    
    # Split the workload across the CPUs.
    pool = multiprocessing.Pool(num_cpus)
    pool_map_async_result = pool.map_async(
            generate_predicted_sedimentation_grid_parallel_pool_function,
            (
                (
                    time,
                    predict_sedimentation_script,
                    scale_sedimentation_rate,
                    mean_age,
                    mean_distance,
                    variance_age,
                    variance_distance,
                    max_age,
                    max_distance,
                    age_distance_polynomial_coefficients,
                    output_dir
                ) for time in range(min_time, max_time + 1, time_step)
            ),
            1) # chunksize

    # Apparently if we use pool.map_async instead of pool.map and then get the results
    # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
    # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    pool_map_async_result.get(99999)
