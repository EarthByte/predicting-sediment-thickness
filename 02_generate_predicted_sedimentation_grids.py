
from ptt.utils.call_system_command import call_system_command
import multiprocessing
import os, shutil
import sys

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
""" ---------- Part 2 of the prediciting-sediment-thickness workflow ----------
This script creates grids (netcdfs) of sediment thickness and sediment rate through time

Requirements amd Inputs:
    - Python
    - Python scripts (predict_sediment_thickness.py, predict sedimentation_rate.py) - should be in this directory
    - GMT 5 (or later)
    - Files associated with a tectonic model, in particular, the agegrids
    - pygplates

To modify the sediment thickness relationship (e.g. for a new present-day agegrid or sediment thickness grid),
you will need to relcalculate the polynomial coefficients, and enter them into this script.

The polynomial coefficients are calculated in the folder 'python_notebooks_and_input_data_archive', in the jupyter
notebooks 'sediment_thick.ipynb' and 'sediment_rate.ipynb'. The values printed by the last cell of this notebook
can be entered into lines 215-222 (for sedimentation rate) and 278-284 (for sediment thickness) of this script.

Outputs:
    - folder named 'sedimentation_output' (or desired name, if changed), with 
      subfolders of sediment thickness and sediment rate grids through time.

2020-02-25: Added comments, created folders within the script itself
2022-08-26: Update parameters for GlobSed and latest agegrids. Modify dirs to be consistent with pt1
"""

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------------------------------------
# --- Set paths and various parameters
# ------------------------------------------

output_base_dir = '.'

# --- distance grids (from part 1)
distance_grid_dir = '%s/distances_1d' % output_base_dir
distance_base_name = 'mean_distance_1d'

# --- output directory name
sediment_output_sub_dir = 'sedimentation_output'

# --- agegrids
agegrid_dir = '/Users/nickywright/Data/Age/Muller2019-Young2019-Cao2020_Agegrids/Muller2019-Young2019-Cao2020_netCDF'   # change folder name if needed
agegrid_filename = 'Muller2019-Young2019-Cao2020_AgeGrid-'    # everything before 'time'
agegrid_filename_ext = 'nc'   # generally 'nc', but sometimes is 'grd'. Do not include the period

# ------------------------------------------
# --- set times and spacing
grid_spacing = 0.1

min_time = 0
max_time = 250
time_step = 1

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

# check if the base output directory exists. If it doesn't, create it.
if not os.path.exists(output_base_dir + '/' + sediment_output_sub_dir):
    print('%s does not exist, creating now... ' % (output_base_dir + '/' + sediment_output_sub_dir))
    os.mkdir(output_base_dir + '/' + sediment_output_sub_dir)


# ----- 
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
    
    py_cmd='python3'
    if os.environ.get('CONDA_PREFIX') or shutil.which('python3') is None:
        py_cmd = 'python'
    
    command_line = [
            py_cmd,
            predict_sedimentation_script,
            '-d',
            '{0}/{1}_{2}.nc'.format(distance_grid_dir, distance_base_name, time),
            '-g',
            '{0}/{1}{2}.{3}'.format(agegrid_dir, agegrid_filename, time, agegrid_filename_ext),
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
    # Only sediment rate requires scaling (sediment thickness does not)...
    if scale_sedimentation_rate is not None:
        command_line.extend([
                '-s',
                str(scale_sedimentation_rate)])
    command_line.extend([
            '--',
            '{0}/sed_{1}_{2}'.format(output_dir, grid_spacing, time)])
    
    #print('Time:', time)
    #print(command_line)

    #print(' '.join(command_line))
    
    # Execute the command.
    call_system_command(command_line)
    
    # Rename the average sedimentation rate and sediment thicknesses files so that 'time' is at the
    # end of the base filename - this way we can import them as time-dependent raster into GPlates.
    for ext in ('xy', 'nc'):
        
        src_sed_rate = '{0}/sed_{1}_{2}_sed_rate.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_sed_rate = '{0}/sed_rate_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        if os.access(dst_sed_rate, os.R_OK):
            os.remove(dst_sed_rate)
        if os.path.exists(src_sed_rate):
            os.rename(src_sed_rate, dst_sed_rate)
        
        src_sed_thick = '{0}/sed_{1}_{2}_sed_thick.{3}'.format(output_dir, grid_spacing, time, ext)
        dst_sed_thick = '{0}/sed_thick_{1}d_{2}.{3}'.format(output_dir, grid_spacing, time, ext)
        
        if os.access(dst_sed_thick, os.R_OK):
            os.remove(dst_sed_thick)
        if os.path.exists(src_sed_thick):
            os.rename(src_sed_thick, dst_sed_thick)


# Wraps around 'generate_predicted_sedimentation_grid()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def generate_predicted_sedimentation_grid_parallel_pool_function(args):
    try:
        return generate_predicted_sedimentation_grid(*args)
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

    times = range(min_time, max_time + 1, time_step)

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

    # updated for GlobSed and TRUNK agegrids
    predict_sedimentation_script = 'predict_sedimentation_rate.py'
    #scale_sedimentation_rate = 1.0  # Keep predicted rate in (cm/Ky).
    scale_sedimentation_rate = 10.0  # Scale predicted rate (cm/Ky) to (m/My).
    mean_age = 61.17716597
    mean_distance = 1835.10750592
    variance_age =  1934.78513885
    variance_distance = 1207587.8548734
    max_age = 191.87276
    max_distance = 3000.
    age_distance_polynomial_coefficients = [
            1.350082937086441, -0.26385415, -0.07516542,  0.39197707, -0.15475392,
        0.        , -0.13196083,  0.02481208, -0.        , -0.47570021]
    
    output_dir = output_base_dir + '/' + sediment_output_sub_dir + '/predicted_rate'
    
    # check if the output dir exists. If not, create
    if not os.path.exists(output_dir):
        print('%s does not exist, creating now... ' % output_dir)
        os.mkdir(output_dir)

    print('Generating predicted sedimentation rate grids...')

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
                    ) for time in times
                ),
                1) # chunksize

        # Apparently if we use pool.map_async instead of pool.map and then get the results
        # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
        # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        pool_map_async_result.get(99999)

    else:
        for time in times:
            generate_predicted_sedimentation_grid(
                    time,
                    predict_sedimentation_script, scale_sedimentation_rate,
                    mean_age, mean_distance,
                    variance_age, variance_distance,
                    max_age, max_distance,
                    age_distance_polynomial_coefficients,
                    output_dir)
    
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

    # updated for GlobSed and TRUNK agegrids (NW 20220826)
    predict_sedimentation_script = 'predict_sediment_thickness.py'
    scale_sedimentation_rate = None  # No scaling - we're predicting sediment thickness (not rate).
    mean_age =  61.18406823
    mean_distance = 1835.28118479
    variance_age = 1934.6999014
    variance_distance = 1207521.8995806
    max_age = 191.87276
    max_distance = 3000.
    age_distance_polynomial_coefficients = [
            5.441401190368497,  0.46893096, -0.07320928, -0.24077496, -0.10840657,
        0.00381672,  0.06831728,  0.01179914,  0.01158149, -0.39880562]
    
    output_dir = output_base_dir + '/' + sediment_output_sub_dir + '/predicted_thickness'

    # check if the output dir exists. If not, create
    if not os.path.exists(output_dir):
        print('%s does not exist, creating now... ' % output_dir)
        os.mkdir(output_dir)
    print('Generating predicted sediment thickness grids...')

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
                    ) for time in times
                ),
                1) # chunksize

        # Apparently if we use pool.map_async instead of pool.map and then get the results
        # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
        # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        pool_map_async_result.get(99999)

    else:
        for time in times:
            generate_predicted_sedimentation_grid(
                    time,
                    predict_sedimentation_script, scale_sedimentation_rate,
                    mean_age, mean_distance,
                    variance_age, variance_distance,
                    max_age, max_distance,
                    age_distance_polynomial_coefficients,
                    output_dir)
