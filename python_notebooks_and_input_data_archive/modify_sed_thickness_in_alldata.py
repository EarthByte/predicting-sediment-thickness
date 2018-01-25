from __future__ import print_function
from call_system_command import call_system_command
import math
import sys

##############################################################################################
# Script to modify the column in 'alldata' file containing mean distance to passive margins. #
##############################################################################################

alldata_input_filename = 'alldata_orig'
alldata_output_filename = 'alldata'

mean_distance_grid_base_filename = 'mean_distance_1.0d_0'
mean_distance_grid_input_filename = mean_distance_grid_base_filename + '.grd'
mean_distance_grid_output_filename = mean_distance_grid_base_filename + '_-180_180_-70_80.grd'


def read_alldata_file(alldata_filename):
    
    alldata = []
    with open(alldata_filename, 'r') as alldata_file:
        for line_number, line in enumerate(alldata_file):

            # Make line number 1-based instead of 0-based.
            line_number = line_number + 1

            # Split the line into strings (separated by whitespace).
            line_string_list = line.split()

            # Need at least two strings per line (for latitude and longitude).
            if len(line_string_list) != 19:
                raise ValueError('Line {0}: Line does not have exactly 19 white-space separated strings.'.format(line_number))

            # Attempt to convert each string into a floating-point number.
            alldata.append(
                    [float(line_string) for line_string in line_string_list])
    
    return alldata


def write_alldata_file(alldata_filename, alldata):
    with open(alldata_filename, 'w') as alldata_file:
        for alldata_line in alldata:
            alldata_file.write(' '.join(str(item) for item in alldata_line) + '\n')


# Restrict latitude range to -70/80 because that corresponds to the extent of the sediment thickness grid data.
# Also using "-T" to convert from pixel to grid registration since distance grids are in pixel registration
# (grid registration gives us pixel values at 1 degree integer lon/lat locations used by other data in 'alldata').
call_system_command(["gmt", "grdsample", "-I1", "-fg", "-T", "-R-180/180/-70/80", mean_distance_grid_input_filename, "-G{0}".format(mean_distance_grid_output_filename)])

# Convert grd to xyz.
mean_distance_output = call_system_command(
        ["gmt", "grd2xyz", mean_distance_grid_output_filename],
        return_stdout=True)

mean_distance_data = []
for line in mean_distance_output.strip().split('\n'):
    line = line.split()
    mean_distance_data.append(
            (float(line[0]), float(line[1]), float(line[2])))
#print(mean_distance_data)

# Create a dictionary mapping lon/lat locations to mean distance.
mean_distance_dict = dict(
        ((float(lon), float(lat)), float(distance))
        for lon, lat, distance in mean_distance_data)


input_alldata = read_alldata_file(alldata_input_filename)

output_alldata = []
for input_data in input_alldata:
    lon, lat = input_data[0], input_data[1]
    distance = mean_distance_dict[(lon, lat)]
    
    # Replace the mean distance to passive margin with the value from mean distance grid.
    output_data = input_data
    output_data[5] = distance
    output_alldata.append(output_data)
    

write_alldata_file(alldata_output_filename, output_alldata)
