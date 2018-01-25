
"""
    Copyright (C) 2017 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


##############################################################################################
# Predict the compacted sediment thickness (metres) at paleo ocean basin point locations and #
# determine average sedimentation rate.                                                      #
##############################################################################################


from __future__ import print_function
import argparse
from call_system_command import call_system_command
import math
import numpy as np
import sys


# Default grid spacing (in degrees) when generating uniform lon/lat spacing of ocean basin points.
DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES = 1

# Default lithology parameters (defaults to shale).
# Obtained from: Bahr, David B., et al. "Exponential approximations to compacted sediment porosity profiles." Computers & Geosciences 27.6 (2001): 691-700.
DEFAULT_SURFACE_POROSITY = 0.63
DEFAULT_POROSITY_EXP_DECAY = 5.71e-4


# Reads the input xy file and returns a list of (lon, lat) points.
def read_input_points(input_points_filename):
    
    input_points = []
    with open(input_points_filename, 'r') as input_points_file:
        for line_number, line in enumerate(input_points_file):

            # Make line number 1-based instead of 0-based.
            line_number = line_number + 1

            # Split the line into strings (separated by whitespace).
            line_string_list = line.split()

            # Need at least two strings per line (for latitude and longitude).
            if len(line_string_list) < 2:
                print('Line {0}: Ignoring point - line does not have at least two white-space separated strings.'.format(
                        line_number), file=sys.stderr)
                continue

            # Attempt to convert each string into a floating-point number.
            try:
                # Use GMT (lon/lat) order.
                lon = float(line_string_list[0])
                lat = float(line_string_list[1])
            except ValueError:
                print('Line {0}: Ignoring point - cannot read lon/lat values.'.format(line_number), file=sys.stderr)
                continue

            input_points.append((lon, lat))
    
    return input_points


def generate_input_points_grid(grid_spacing_degrees):
    
    if grid_spacing_degrees == 0:
        return
    
    input_points = []
    
    # Data points start *on* dateline (-180).
    # If 180 is an integer multiple of grid spacing then final longitude also lands on dateline (+180).
    num_latitudes = int(math.floor(180.0 / grid_spacing_degrees)) + 1
    num_longitudes = int(math.floor(360.0 / grid_spacing_degrees)) + 1
    for lat_index in range(num_latitudes):
        lat = -90 + lat_index * grid_spacing_degrees
        
        for lon_index in range(num_longitudes):
            lon = -180 + lon_index * grid_spacing_degrees
            
            input_points.append((lon, lat))
    
    return (input_points, num_longitudes, num_latitudes)


# Returns a list of scalars (one per (lon, lat) point in the 'input_points' list).
# For input points outside the scalar grid then scalars will be Nan (ie, 'math.isnan(scalar)' will return True).
def get_positions_and_scalars(input_points, scalar_grid_filename, max_scalar=None):
    
    input_points_data = ''.join('{0} {1}\n'.format(lon, lat) for lon, lat in input_points)

    # The command-line strings to execute GMT 'grdtrack'.
    grdtrack_command_line = ["gmt", "grdtrack", "-nl", "-G{0}".format(scalar_grid_filename)]
    stdout_data = call_system_command(grdtrack_command_line, stdin=input_points_data, return_stdout=True)
    
    lon_lat_scalar_list = []
    
    # Read lon, lat and scalar values from the output of 'grdtrack'.
    for line in stdout_data.splitlines():
        if line.strip().startswith('#'):
            continue
        
        line_data = line.split()
        num_values = len(line_data)
        
        # If just a line containing white-space then skip to next line.
        if num_values == 0:
            continue
        
        if num_values < 3:
            print('Ignoring line "{0}" - has fewer than 3 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            lon = float(line_data[0])
            lat = float(line_data[1])
            
            # The scalar got appended to the last column by 'grdtrack'.
            scalar = float(line_data[-1])
            
            # If the point is outside the grid then the scalar grid will return 'NaN'.
            if math.isnan(scalar):
                #print('Ignoring line "{0}" - point is outside scalar grid.'.format(line), file=sys.stderr)
                continue
            
            # Clamp to max value if requested.
            if (max_scalar is not None and
                scalar > max_scalar):
                scalar = max_scalar
            
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point lon, lat and scalar values.'.format(line), file=sys.stderr)
            continue
        
        lon_lat_scalar_list.append((lon, lat, scalar))
    
    return lon_lat_scalar_list


def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def write_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, num_grid_longitudes, num_grid_latitudes):
    
    # The command-line strings to execute GMT 'nearneighbor'.
    # For example "nearneighbor output_mean_distance.xy -R-180/180/-90/90 -I1 -N4 -S1d -Goutput_mean_distance.grd=cf"
    # ...with "=cf" generating NetCDF-3, instead of NetCDF-4, files (GPlates 2.0 can only load NetCDF-3).
    gmt_command_line = [
            "gmt",
            "nearneighbor",
            xyz_filename.encode(sys.getfilesystemencoding()),
            "-N4",
            "-S{0}d".format(1.5 * grid_spacing),
            "-I{0}".format(grid_spacing),
            "-R{0}/{1}/{2}/{3}".format(-180, 180, -90, 90),
            # Use GMT gridline registration since our input point grid has data points on the grid lines.
            # Gridline registration is the default so we don't need to force pixel registration...
            #"-r", # Force pixel registration since data points are at centre of cells.
            "-G{0}=cf".format(grd_filename.encode(sys.getfilesystemencoding()))]
    call_system_command(gmt_command_line)


# Decompact total sediment thickness assuming a single lithology of the specified surface porosity and exponential decay.
#
# The decompacted (ie, surface) thickness d(z) of a compacted layer of height t(z) at depth z is:
#
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# ...due to the fact that the grain volume doesn't change with depth:
#
#   (1 - porosity(z1)) * t(z1) = (1 - porosity(z2)) * t(z2)
#
# ...and substituting z1=z and z2=0, and using d(z) = t(0) gives:
#
#   (1 - porosity(z)) * t(z) = (1 - porosity(0)) * d(z)
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# Where the porosity is:
#
#   porosity(z) = porosity(0) * exp(-c * z)
#
# The total decompacted thickness is:
#
#   D = Integral(d(z), z = 0 -> T)
#
# ...where T is the total compacted thickness.
#
#   D = Integral(dz * (1 - porosity(z)) / (1 - porosity(0)), z = 0 -> T)
#   D = Integral(dz * (1 - porosity(0) * exp(-c * z)) / (1 - porosity(0)), z = 0 -> T)
#
#     = (T + (porosity(0) / c) * (exp(-c * T) - 1)) / (1 - porosity(0))
#
def decompact_sediment_thickness(thickness, surface_porosity, porosity_exp_decay):
    return (thickness + (surface_porosity / porosity_exp_decay) * (math.exp(-porosity_exp_decay * thickness) - 1)) / (1 - surface_porosity)


def predict_sediment_thickness(
        age,
        distance,
        mean_age,
        mean_distance,
        std_deviation_age,
        std_deviation_distance,
        age_distance_polynomial_coefficients):
    
    # We need to remove the mean and scale to unit variance (for the age and distance values)
    # based on the machine learning training scaler.
    # See http://scikit-learn.org/stable/modules/preprocessing.html#standardization-or-mean-removal-and-variance-scaling
    age = (age - mean_age) / std_deviation_age
    distance = (distance - mean_distance) / std_deviation_distance
    
    age_distance_polynomial_features = [
            1, age, distance, age*age, age*distance, distance*distance, age*age*age, age*age*distance, age*distance*distance, distance*distance*distance]
    
    # Evaluate the polynomial to get the log of the predicated sediment thickness.
    log_sediment_thickness = sum(age_distance_polynomial_coefficients[i] * age_distance_polynomial_features[i] for i in range(10))
    
    # Return the predicted sediment thickness (not as a logarithm).
    return math.exp(log_sediment_thickness)


def predict_sedimentation(
        input_points, # List of (lon, lat) tuples,
        age_grid_filename,
        distance_grid_filename,
        mean_age,
        mean_distance,
        variance_age,
        variance_distance,
        age_distance_polynomial_coefficients,
        surface_porosity,
        porosity_exp_decay,
        max_age = None,
        max_distance = None,
        sedimentation_rate_scale = 1.0):
    """
    Predicts compacted sediment thickness for each ocean basin point based on its age and its
    mean distance to passive continental margins. Also determines average sedimentation rate.
    
    Returns: A 2-tuple of (lon_lat_average_sedimentation_rate_list, lon_lat_sediment_thickness_list)
             where the two lists contain 3-tuples of (lon, lat, value).
    """
    
    # Get the input point ages and mean distances to passive continental margins.
    lon_lat_age_list = get_positions_and_scalars(input_points, age_grid_filename, max_age)
    lon_lat_distance_list = get_positions_and_scalars(input_points, distance_grid_filename, max_distance)
    if not lon_lat_age_list or not lon_lat_distance_list:
        return
    
    # Merge the age and distance lists.
    # Only keep points where there are age *and* distance values.
    lon_lat_age_distance_list = []
    age_dict = dict(((lon, lat), age) for lon, lat, age in lon_lat_age_list)
    for lon, lat, distance in lon_lat_distance_list:
        if (lon, lat) in age_dict:
            age = age_dict[(lon, lat)]
            lon_lat_age_distance_list.append((lon, lat, age, distance))
    
    #
    # Calculate mean/variance statistics on the age/distance data.
    #
    # Note: We don't use the mean/variance of data (instead need to use mean/variance from trained data).
    #
    #sum_ages = 0
    #sum_square_ages = 0
    #sum_distances = 0
    #sum_square_distances = 0
    #for lon, lat, age, distance in lon_lat_age_distance_list:
    #    sum_ages += age
    #    sum_square_ages += age * age
    #    sum_distances += distance
    #    sum_square_distances += distance * distance
    #num_age_distance_points = len(lon_lat_age_distance_list)
    #mean_age = sum_ages / num_age_distance_points
    #mean_distance = sum_distances / num_age_distance_points
    #variance_age = (sum_square_ages / num_age_distance_points) - (mean_age * mean_age)
    #variance_distance = (sum_square_distances / num_age_distance_points) - (mean_distance * mean_distance)
    
    std_deviation_age = math.sqrt(variance_age)
    std_deviation_distance = math.sqrt(variance_distance)
    
    # For each ocean basin point predict compacted sediment thickness and determine the average sedimentation rate.
    lon_lat_average_sedimentation_rate_list = []
    lon_lat_sediment_thickness_list = []
    for lon, lat, age, distance in lon_lat_age_distance_list:
        
        predicted_sediment_thickness = predict_sediment_thickness(
                age, distance,
                mean_age, mean_distance,
                std_deviation_age, std_deviation_distance,
                age_distance_polynomial_coefficients)
        lon_lat_sediment_thickness_list.append((lon, lat, predicted_sediment_thickness))
        
        # Decompact the sediment thickness and divide by ocean floor age.
        #
        # Adhoc dealing of problems near mid-ocean ridges where age is very young.
        # Avoids divide-by-zero and very large sedimentation rates.
        if age < 1:
            age = 1
        average_sedimentation_rate = sedimentation_rate_scale * decompact_sediment_thickness(
                predicted_sediment_thickness, surface_porosity, porosity_exp_decay) / age
        lon_lat_average_sedimentation_rate_list.append((lon, lat, average_sedimentation_rate))
    
    return lon_lat_average_sedimentation_rate_list, lon_lat_sediment_thickness_list
    
    
def write_sediment_data(
        average_sedimentation_rate_data,
        sediment_thickness_data,
        output_filename_prefix,
        output_filename_extension,
        output_grd_file = None):
    
    average_sedimentation_rate_suffix = 'sed_rate'
    sediment_thickness_suffix = 'sed_thick'
    
    # Write average sedimentation rate and sediment thickness XYZ files.
    average_sedimentation_rate_xyz_filename = u'{0}_{1}.{2}'.format(
            output_filename_prefix, average_sedimentation_rate_suffix, output_filename_extension)
    sediment_thickness_xyz_filename = u'{0}_{1}.{2}'.format(
            output_filename_prefix, sediment_thickness_suffix, output_filename_extension)
    
    write_xyz_file(average_sedimentation_rate_xyz_filename, average_sedimentation_rate_data)
    write_xyz_file(sediment_thickness_xyz_filename, sediment_thickness_data)
    
    # Write average sedimentation rate and sediment thickness GRD files.
    if output_grd_file:
        grid_spacing, num_grid_longitudes, num_grid_latitudes = output_grd_file
        
        average_sedimentation_rate_grd_filename = u'{0}_{1}.grd'.format(
                output_filename_prefix, average_sedimentation_rate_suffix)
        sediment_thickness_grd_filename = u'{0}_{1}.grd'.format(
                output_filename_prefix, sediment_thickness_suffix)
        
        write_grd_file_from_xyz(
                average_sedimentation_rate_grd_filename,
                average_sedimentation_rate_xyz_filename,
                grid_spacing, num_grid_longitudes, num_grid_latitudes)
        write_grd_file_from_xyz(
                sediment_thickness_grd_filename,
                sediment_thickness_xyz_filename,
                grid_spacing, num_grid_longitudes, num_grid_latitudes)


if __name__ == '__main__':
    
    __description__ = \
    """Find the average sedimentation rate (metres / million years) at ocean basin point locations assuming a single lithology for the entire sediment thickness.
    
    The sediment thicknesses are decompacted and then divided by age at ocean basin point locations to get average sedimentation rates.
    
    An input xy file can be specified. If one is not specified then a uniform lon/lat grid of points will be generated internally
    (its grid spacing is controlled with the '-i' option and the grid point latitudes/longitudes use GMT gridline registration,
    eg, for a 1 degree spacing the latitudes are -90, -89, ..., 90). If an input xy file is specified then it can contain
    arbitrary (lon, lat) point locations (it can also contain extra 3rd, 4th, etc, columns but they are ignored).
    
    By default a single 'xy' file containing average sedimentation rates is output.
    An optional 'grd' can also be outpu (see the '-w' option)
    
    If an ocean basin point falls outside the age grid or the sediment thickness grid then it is ignored.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -d distance.grd -g age_grid.nc -i 1 -w -- sediment
     """
    
    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-d', '--distance_grid_filename', type=str, required=True,
            metavar='distance_grid_filename',
            help='The distance grid filename containing the mean distances, to passive continental margins, of each '
                'ocean basin point (reconstructed from its time of appearance to the paleo time of the age grid).')
    parser.add_argument('-g', '--age_grid_filename', type=str, required=True,
            metavar='age_grid_filename',
            help='The age grid filename used to find the age of each ocean basin point.')
    parser.add_argument('-w', '--output_grd_file', action='store_true',
            help='Also generate a grd file. '
                'By default only an xyz file is written. '
                'Can only be specified if "ocean_basin_points_filename" is not specified '
                '(ie, ocean basin points must be on a uniform lon/lat grid).')
    
    parser.add_argument('-i', '--ocean_basin_grid_spacing', type=float,
            help='The grid spacing (in degrees) of ocean basin points in lon/lat space. '
                'The grid point latitudes/longitudes are offset by half the grid spacing '
                '(eg, for a 1 degree spacing the latitudes are -89.5, -88.5, ..., 89.5). '
                'Can only be specified if "ocean_basin_points_filename" is not specified. '
                'Defaults to {0} degrees.'.format(
                        DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES))
    
    parser.add_argument('-s', '--sedimentation_rate_scale', type=float, default=1.0,
            metavar='sedimentation_rate_scale',
            help='Scale factor to apply to average sedimentation rate to convert it *from* metres/My. '
                'For example, if you want rate in units cm/Ky then set the scale factor to 0.1. '
                'Defaults to no scaling (1.0).')
    
    parser.add_argument('-m', '--mean_age_distance', type=float, nargs=2, required=True,
            metavar=('mean_age', 'mean_distance'),
            help='The mean of age and distance as two consecutive values. '
                'These values come from the machine learning training scaler and are used to '
                'remove the mean (for age and distance).')
    parser.add_argument('-v', '--variance_age_distance', type=float, nargs=2, required=True,
            metavar=('variance_age', 'variance_distance'),
            help='The variance of age and distance as two consecutive values. '
                'These values come from the machine learning training scaler and are used to '
                'scale to unit variance (for age and distance).')
    parser.add_argument('-x', '--clamp_age_distance', type=float, nargs=2, default=None,
            metavar=('max_age', 'max_distance'),
            help='Optional maximum of age and distance as two consecutive values. '
                'These values come from the machine learning training stage and are used to clamp the '
                'age and distance (because values above these are not represented well in trained data). '
                'Defaults to no clamping (of age or distance).')
    parser.add_argument('-f', '--age_distance_polynomial_coefficients', type=float, nargs=10, required=True,
            metavar=('constant', 'age', 'distance', 'age*age', 'age*distance', 'distance*distance',
                    'age*age*age', 'age*age*distance', 'age*distance*distance', 'distance*distance*distance'),
            help='The polynomial coefficients used to predict sedimentation rate from age and distance. '
                'There should be ten consecutive values representing the polynomial features: '
                'constant, age, distance, age*age, age*distance, distance*distance, age*age*age, age*age*distance, '
                'age*distance*distance and distance*distance*distance. '
                'These values come from the machine learning training scaler.')
    
    parser.add_argument('-p', '--surface_porosity', type=float, default=DEFAULT_SURFACE_POROSITY,
            help='The surface porosity in range [0,1]. '
                'Defaults to {0} for shale.'.format(DEFAULT_SURFACE_POROSITY))
    parser.add_argument('-c', '--porosity_exp_decay', type=float, default=DEFAULT_POROSITY_EXP_DECAY,
            help='The exponential decay of porosity with depth where '
            '"porosity = surface_porosity * exp(-porosity_exp_decay * depth)". '
                'Defaults to {0} for shale.'.format(DEFAULT_POROSITY_EXP_DECAY))
    
    def parse_unicode(value_string):
        try:
            # Filename uses the system encoding - decode from 'str' to 'unicode'.
            filename = value_string.decode(sys.getfilesystemencoding())
        except UnicodeDecodeError:
            raise argparse.ArgumentTypeError("Unable to convert filename %s to unicode" % value_string)
        
        return filename
    
    parser.add_argument('ocean_basin_points_filename', type=parse_unicode, nargs='?',
            metavar='ocean_basin_points_filename',
            help='Optional input xy file containing the ocean basin point locations. '
                'If not specified then a uniform lon/lat grid of points is generated. '
                'Can only be specified if "ocean_basin_grid_spacing" and "output_grd_file" are not specified.')
    
    parser.add_argument('output_filename_prefix', type=parse_unicode,
            metavar='output_filename_prefix',
            help='The output xy filename prefix used in all output filenames.')
    parser.add_argument('-e', '--output_filename_extension', type=str, default='xy',
            metavar='output_filename_extension',
            help='The output xy filename extension. Defaults to "xy".')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    if args.ocean_basin_points_filename is not None:
        if args.ocean_basin_grid_spacing is not None:
            raise argparse.ArgumentTypeError(
                "'ocean_basin_grid_spacing' and 'ocean_basin_points_filename' cannot both be specified.")
        if args.output_grd_file is not None:
            raise argparse.ArgumentTypeError(
                "'output_grd_file' and 'ocean_basin_points_filename' cannot both be specified.")
    
    # Get the input points.
    if args.ocean_basin_points_filename is not None:
        input_points = read_input_points(args.ocean_basin_points_filename)
    else:
        if args.ocean_basin_grid_spacing is None:
            args.ocean_basin_grid_spacing = DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES
        input_points, num_grid_longitudes, num_grid_latitudes = generate_input_points_grid(args.ocean_basin_grid_spacing)
    
    if input_points is None:
        sys.exit(1)
    
    average_sedimentation_rate_data, sediment_thickness_data = predict_sedimentation(
            input_points,
            args.age_grid_filename,
            args.distance_grid_filename,
            args.mean_age_distance[0], # mean_age
            args.mean_age_distance[1], # mean_distance
            args.variance_age_distance[0], # variance_age
            args.variance_age_distance[1], # variance_distance
            args.age_distance_polynomial_coefficients,
            args.surface_porosity,
            args.porosity_exp_decay,
            args.clamp_age_distance[0] if args.clamp_age_distance is not None else None, # max_age
            args.clamp_age_distance[1] if args.clamp_age_distance is not None else None, # max_distance
            args.sedimentation_rate_scale)

    if average_sedimentation_rate_data is None:
        sys.exit(1)
    
    write_sediment_data(
            average_sedimentation_rate_data,
            sediment_thickness_data,
            args.output_filename_prefix,
            args.output_filename_extension,
            (args.ocean_basin_grid_spacing, num_grid_longitudes, num_grid_latitudes) if args.output_grd_file else None)
    
    sys.exit(0)
