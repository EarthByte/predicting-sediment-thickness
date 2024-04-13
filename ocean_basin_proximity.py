
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


#############################################################################################################################
# Find the minimum distance of ocean basin point locations to topological boundaries or non-topological features over time. #
#############################################################################################################################



import argparse
import math
import multiprocessing
import numpy as np
import os
# Try importing 'ptt' first. If that fails then try 'gplately.ptt' (GPlately now contains PlateTectonicTools).
try:
    from ptt.utils.call_system_command import call_system_command
    import ptt.utils.proximity_query as proximity_query
except ImportError:
    from gplately.ptt.utils.call_system_command import call_system_command
    import gplately.ptt.utils.proximity_query as proximity_query
import pygplates
from scipy.spatial import KDTree
import shortest_path
import sys
import time as time_profile


# Enable CPU/memory profiling.
ENABLE_CPU_PROFILING = False
ENABLE_MEMORY_PROFILING = False


# Class to profile CPU usage.
class CpuProfile(object):
    def __init__(self, enable_profiling=False):
        self.enable_profiling = enable_profiling
    
    def start_proximity(self):
        """Call at the start of proximity()."""
        if self.enable_profiling:
            self.time_snapshot_start_proximity = self.profile()

            self.time_usage_proximity = 0.0
            self.time_usage_read_input_data = 0.0
            self.time_usage_reconstruct_and_calculate_distances = 0.0
            self.time_usage_read_age_grid = 0.0
            self.time_usage_reconstruct_proximity = 0.0
            self.time_usage_calculate_distances = 0.0
            self.time_usage_obstacle_reconstruct_resolve = 0.0
            self.time_usage_create_obstacle_grids = 0.0
            self.time_usage_calculate_obstacle_distances = 0.0
            self.time_usage_obstacle_create_obstacle_grid = 0.0
            self.time_usage_obstacle_create_distance_grid = 0.0
            self.time_usage_obstacle_calc_distances = 0.0
            self.time_usage_reconstruct_time_step = 0.0
    def end_proximity(self):
        """Call at the end of proximity()."""
        if self.enable_profiling:
            self.time_usage_proximity += self.profile() - self.time_snapshot_start_proximity
    
    def start_read_input_data(self):
        if self.enable_profiling:
            self.time_snapshot_start_read_input_data = self.profile()
    def end_read_input_data(self):
        if self.enable_profiling:
            self.time_usage_read_input_data += self.profile() - self.time_snapshot_start_read_input_data
    
    def start_reconstruct_and_calculate_distances(self):
        if self.enable_profiling:
            self.time_snapshot_reconstruct_and_calculate_distances = self.profile()
    def end_reconstruct_and_calculate_distances(self):
        if self.enable_profiling:
            self.time_usage_reconstruct_and_calculate_distances += self.profile() - self.time_snapshot_reconstruct_and_calculate_distances
    
    def start_read_age_grid(self):
        if self.enable_profiling:
            self.time_snapshot_read_age_grid = self.profile()
    def end_read_age_grid(self):
        if self.enable_profiling:
            self.time_usage_read_age_grid += self.profile() - self.time_snapshot_read_age_grid
    
    def start_reconstruct_proximity(self):
        if self.enable_profiling:
            self.time_snapshot_start_reconstruct_proximity = self.profile()
    def end_reconstruct_proximity(self):
        if self.enable_profiling:
            self.time_usage_reconstruct_proximity += self.profile() - self.time_snapshot_start_reconstruct_proximity
    
    def start_calculate_distances(self):
        if self.enable_profiling:
            self.time_snapshot_start_calculate_distances = self.profile()
    def end_calculate_distances(self):
        if self.enable_profiling:
            self.time_usage_calculate_distances += self.profile() - self.time_snapshot_start_calculate_distances
    
    def start_obstacle_reconstruct_resolve(self):
        if self.enable_profiling:
            self.time_snapshot_start_obstacle_reconstruct_resolve = self.profile()
    def end_obstacle_reconstruct_resolve(self):
        if self.enable_profiling:
            self.time_usage_obstacle_reconstruct_resolve += self.profile() - self.time_snapshot_start_obstacle_reconstruct_resolve
    
    def start_create_obstacle_grids(self):
        if self.enable_profiling:
            self.time_snapshot_start_create_obstacle_grids = self.profile()
    def end_create_obstacle_grids(self):
        if self.enable_profiling:
            self.time_usage_create_obstacle_grids += self.profile() - self.time_snapshot_start_create_obstacle_grids
    
    def start_calculate_obstacle_distances(self):
        if self.enable_profiling:
            self.time_snapshot_start_calculate_obstacle_distances = self.profile()
    def end_calculate_obstacle_distances(self):
        if self.enable_profiling:
            self.time_usage_calculate_obstacle_distances += self.profile() - self.time_snapshot_start_calculate_obstacle_distances
    
    def start_reconstruct_time_step(self):
        if self.enable_profiling:
            self.time_snapshot_start_reconstruct_time_step = self.profile()
    def end_reconstruct_time_step(self):
        if self.enable_profiling:
            self.time_usage_reconstruct_time_step += self.profile() - self.time_snapshot_start_reconstruct_time_step
    

    def start_write_proximity_data(self):
        """Call at the start of write_proximity_data()."""
        if self.enable_profiling:
            self.time_snapshot_start_write_proximity_data = self.profile()

            self.time_usage_write_proximity_data = 0.0
            self.time_usage_write_xyz_file = 0.0
            self.time_usage_write_grd_file_from_xyz = 0.0
            self.time_usage_upscaled_mask_generate_input_points = 0.0
            self.time_usage_calculate_upscaled_mask_interpolation_params = 0.0
            self.time_usage_upscaled_sample_age_grid = 0.0
            self.time_usage_extract_upscale_mask_samples = 0.0
            self.time_usage_create_and_query_kdtrees = 0.0
            self.time_usage_calc_interpolation_weights = 0.0
            self.time_usage_write_upscaled_grd_file = 0.0
            self.time_usage_calc_upscaled_scalars = 0.0
    def end_write_proximity_data(self):
        """Call at the end of write_proximity_data()."""
        if self.enable_profiling:
            self.time_usage_write_proximity_data += self.profile() - self.time_snapshot_start_write_proximity_data
    
    def start_write_xyz_file(self):
        if self.enable_profiling:
            self.time_snapshot_start_write_xyz_file = self.profile()
    def end_write_xyz_file(self):
        if self.enable_profiling:
            self.time_usage_write_xyz_file += self.profile() - self.time_snapshot_start_write_xyz_file
    
    def start_write_grd_file_from_xyz(self):
        if self.enable_profiling:
            self.time_snapshot_start_write_grd_file_from_xyz = self.profile()
    def end_write_grd_file_from_xyz(self):
        if self.enable_profiling:
            self.time_usage_write_grd_file_from_xyz += self.profile() - self.time_snapshot_start_write_grd_file_from_xyz
    
    def start_upscaled_mask_generate_input_points(self):
        if self.enable_profiling:
            self.time_snapshot_start_upscaled_mask_generate_input_points = self.profile()
    def end_upscaled_mask_generate_input_points(self):
        if self.enable_profiling:
            self.time_usage_upscaled_mask_generate_input_points += self.profile() - self.time_snapshot_start_upscaled_mask_generate_input_points
    
    def start_calculate_upscaled_mask_interpolation_params(self):
        if self.enable_profiling:
            self.time_snapshot_start_calculate_upscaled_mask_interpolation_params = self.profile()
    def end_calculate_upscaled_mask_interpolation_params(self):
        if self.enable_profiling:
            self.time_usage_calculate_upscaled_mask_interpolation_params += self.profile() - self.time_snapshot_start_calculate_upscaled_mask_interpolation_params
    
    def start_upscaled_sample_age_grid(self):
        if self.enable_profiling:
            self.time_snapshot_start_upscaled_sample_age_grid = self.profile()
    def end_upscaled_sample_age_grid(self):
        if self.enable_profiling:
            self.time_usage_upscaled_sample_age_grid += self.profile() - self.time_snapshot_start_upscaled_sample_age_grid
    
    def start_extract_upscale_mask_samples(self):
        if self.enable_profiling:
            self.time_snapshot_start_extract_upscale_mask_samples = self.profile()
    def end_extract_upscale_mask_samples(self):
        if self.enable_profiling:
            self.time_usage_extract_upscale_mask_samples += self.profile() - self.time_snapshot_start_extract_upscale_mask_samples
    
    def start_create_and_query_kdtrees(self):
        if self.enable_profiling:
            self.time_snapshot_start_create_and_query_kdtrees = self.profile()
    def end_create_and_query_kdtrees(self):
        if self.enable_profiling:
            self.time_usage_create_and_query_kdtrees += self.profile() - self.time_snapshot_start_create_and_query_kdtrees
    
    def start_calc_interpolation_weights(self):
        if self.enable_profiling:
            self.time_snapshot_start_calc_interpolation_weights = self.profile()
    def end_calc_interpolation_weights(self):
        if self.enable_profiling:
            self.time_usage_calc_interpolation_weights += self.profile() - self.time_snapshot_start_calc_interpolation_weights
    
    def start_write_upscaled_grd_file(self):
        if self.enable_profiling:
            self.time_snapshot_start_write_upscaled_grd_file = self.profile()
    def end_write_upscaled_grd_file(self):
        if self.enable_profiling:
            self.time_usage_write_upscaled_grd_file += self.profile() - self.time_snapshot_start_write_upscaled_grd_file
    
    def start_calc_upscaled_scalars(self):
        if self.enable_profiling:
            self.time_snapshot_start_calc_upscaled_scalars = self.profile()
    def end_calc_upscaled_scalars(self):
        if self.enable_profiling:
            self.time_usage_calc_upscaled_scalars += self.profile() - self.time_snapshot_start_calc_upscaled_scalars
    

    def print_usage(self, age_grid_paleo_times):
        """Call to print CPU usage."""
        if self.enable_profiling:
            scale_to_seconds = 1e-9  # convert nanoseconds to seconds
            print( "proximity() CPU usage:")
            print(f"  Age grid paleo times: {age_grid_paleo_times}")
            print(f"    Proximity: {self.time_usage_proximity * scale_to_seconds:.2f} seconds")
            print(f"      Read input data: {self.time_usage_read_input_data * scale_to_seconds:.2f} seconds")
            print(f"      Reconstruct and calculate distances: {self.time_usage_reconstruct_and_calculate_distances * scale_to_seconds:.2f} seconds")
            print(f"        Read age grid: {self.time_usage_read_age_grid * scale_to_seconds:.2f} seconds")
            print(f"        Reconstruct proximity: {self.time_usage_reconstruct_proximity * scale_to_seconds:.2f} seconds")
            print(f"        Calculate distances: {self.time_usage_calculate_distances * scale_to_seconds:.2f} seconds")
            print(f"          Obstacle reconstruct/resolve: {self.time_usage_obstacle_reconstruct_resolve * scale_to_seconds:.2f} seconds")
            print(f"          Obstacle create obstacle grids: {self.time_usage_create_obstacle_grids * scale_to_seconds:.2f} seconds")
            print(f"          Obstacle calculate distances: {self.time_usage_calculate_obstacle_distances * scale_to_seconds:.2f} seconds")
            print(f"        Reconstruct time steps: {self.time_usage_reconstruct_time_step * scale_to_seconds:.2f} seconds")
            print(f"    WriteProximityData: {self.time_usage_write_proximity_data * scale_to_seconds:.2f} seconds")
            print(f"      Write xyz file: {self.time_usage_write_xyz_file * scale_to_seconds:.2f} seconds")
            print(f"      Write grd file from xyz: {self.time_usage_write_grd_file_from_xyz * scale_to_seconds:.2f} seconds")
            print(f"      Generate upscaled mask input points: {self.time_usage_upscaled_mask_generate_input_points * scale_to_seconds:.2f} seconds")
            print(f"      Get upscaled mask interpolation parameters: {self.time_usage_calculate_upscaled_mask_interpolation_params * scale_to_seconds:.2f} seconds")
            print(f"        Sample age grid at upscaled grid spacing: {self.time_usage_upscaled_sample_age_grid * scale_to_seconds:.2f} seconds")
            print(f"        Extract upscaled mask samples: {self.time_usage_extract_upscale_mask_samples * scale_to_seconds:.2f} seconds")
            print(f"        Create and query k-d trees: {self.time_usage_create_and_query_kdtrees * scale_to_seconds:.2f} seconds")
            print(f"        Calculate interpolation weights: {self.time_usage_calc_interpolation_weights * scale_to_seconds:.2f} seconds")
            print(f"      Write upscaled grd file: {self.time_usage_write_upscaled_grd_file * scale_to_seconds:.2f} seconds")
            print(f"        Calculate upscaled scalars: {self.time_usage_calc_upscaled_scalars * scale_to_seconds:.2f} seconds")
    
    @staticmethod
    def profile():
        return time_profile.perf_counter_ns()  # in nanoseconds

# Profile CPU usage (currently just the profile() function).
cpu_profile = CpuProfile(ENABLE_CPU_PROFILING)


from itertools import chain
from collections import deque
try:
    from reprlib import repr
except ImportError:
    pass

# Class to profile memory usage.
class MemoryProfile(object):
    def __init__(self, enable_profiling=False):
        self.enable_profiling = enable_profiling
    
    def print_object_memory_usage(self, obj, obj_name, decimal_places=2):
        """Print the total memory usage of an object (in MB)."""
        if self.enable_profiling:
            obj_memory_usage = self.total_memory_size(obj) / 1e6  # in MB
            print('Memory usage "{}": {:.{}f}MB'.format(obj_name, obj_memory_usage, decimal_places))

    # This function was obtained from https://code.activestate.com/recipes/577504/
    @staticmethod
    def total_memory_size(o, handlers={}, verbose=False):
        """ Returns the approximate memory footprint an object and all of its contents.

        Automatically finds the contents of the following builtin containers and
        their subclasses:  tuple, list, deque, dict, set and frozenset.
        To search other containers, add handlers to iterate over their contents:

            handlers = {SomeContainerClass: iter,
                        OtherContainerClass: OtherContainerClass.get_elements}

        """
        dict_handler = lambda d: chain.from_iterable(d.items())
        all_handlers = {tuple: iter,
                        list: iter,
                        deque: iter,
                        dict: dict_handler,
                        set: iter,
                        frozenset: iter,
                    }
        all_handlers.update(handlers)     # user handlers take precedence
        seen = set()                      # track which object id's have already been seen
        default_size = sys.getsizeof(0)       # estimate sizeof object without __sizeof__

        def sizeof(o):
            if id(o) in seen:       # do not double count the same object
                return 0
            seen.add(id(o))
            s = sys.getsizeof(o, default_size)

            if verbose:
                print(s, type(o), repr(o), file=sys.stderr)

            for typ, handler in all_handlers.items():
                if isinstance(o, typ):
                    s += sum(map(sizeof, handler(o)))
                    break
            else:
                if not hasattr(o.__class__, '__slots__'):
                    if hasattr(o, '__dict__'):
                        s+=sizeof(o.__dict__) # no __slots__ *usually* means a __dict__, but some special builtin classes (such as `type(None)`) have neither
                    # else, `o` has no attributes at all, so sys.getsizeof() actually returned the correct value
                else:
                    s+=sum(sizeof(getattr(o, x)) for x in o.__class__.__slots__ if hasattr(o, x))
            return s

        return sizeof(o)

# Profile memory usage (currently just the profile() function).
memory_profile = MemoryProfile(ENABLE_MEMORY_PROFILING)


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
                print('WARNING: Line {}: Ignoring point - line does not have at least two white-space separated strings.'.format(
                        line_number), file=sys.stderr)
                continue

            # Attempt to convert each string into a floating-point number.
            try:
                # Use GMT (lon/lat) order.
                lon = float(line_string_list[0])
                lat = float(line_string_list[1])
            except ValueError:
                print('WARNING: Line {}: Ignoring point - cannot read lon/lat values.'.format(line_number), file=sys.stderr)
                continue

            input_points.append((lon, lat))
    
    return np.array(input_points)  # numpy array uses less memory


def generate_input_points_grid(grid_spacing_degrees):
    
    if grid_spacing_degrees == 0:
        raise ValueError('Grid spacing cannot be zero.')
    
    # Data points start *on* dateline (-180).
    # If 180 is an integer multiple of grid spacing then final longitude also lands on dateline (+180).
    num_latitudes = int(math.floor(180.0 / grid_spacing_degrees)) + 1
    num_longitudes = int(math.floor(360.0 / grid_spacing_degrees)) + 1

    # Generate the input points on the grid.
    input_points_mesh_grid = np.meshgrid(
            np.linspace(-180, 180, num_longitudes),
            np.linspace(-90, 90, num_latitudes))
    input_points = np.array(input_points_mesh_grid).reshape(2, -1).T
    
    return (input_points, num_longitudes, num_latitudes)


# Returns a list of ages (one per (lon, lat) point in the 'input_points' list).
# Input points outside the age grid are ignored.
def get_positions_and_ages(input_points, age_grid_filename):
    
    input_points_data = ''.join('{} {}\n'.format(lon, lat) for lon, lat in input_points)
    
    stdout_data = call_system_command(
            # The command-line strings to execute GMT 'grdtrack'...
            ["gmt", "grdtrack", "-G{}".format(age_grid_filename)],
            stdin=input_points_data,
            return_stdout=True)
    
    #print('Stdout: {}'.format(stdout_data))
    
    lon_lat_age_list = []
    
    # Read lon, lat and age values from the output of 'grdtrack'.
    for line in stdout_data.splitlines():
        if line.strip().startswith('#'):
            continue
        
        line_data = line.split()
        num_values = len(line_data)
        
        # If just a line containing white-space then skip to next line.
        if num_values == 0:
            continue
        
        if num_values < 3:
            print('WARNING: Ignoring line "{}" - has fewer than 3 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            lon = float(line_data[0])
            lat = float(line_data[1])
            
            # The age got appended to the last column by 'grdtrack'.
            age = float(line_data[-1])
            
            # If the point is outside the ocean basin region then the age grid will return 'NaN'.
            if math.isnan(age):
                #print('WARNING: Ignoring line "{}" - point is outside ocean basin (age grid).'.format(line), file=sys.stderr)
                continue
            
        except ValueError:
            print('WARNING: Ignoring line "{}" - cannot read floating-point lon, lat and age values.'.format(line), file=sys.stderr)
            continue
        
        lon_lat_age_list.append((lon, lat, age))
    
    return lon_lat_age_list


def write_xyz_file(output_filename, output_data):
    cpu_profile.start_write_xyz_file()

    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')
    
    cpu_profile.end_write_xyz_file()


def write_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, use_nearneighbor = True):
    cpu_profile.start_write_grd_file_from_xyz()
    
    if use_nearneighbor:
        # The command-line strings to execute GMT 'nearneighbor'.
        # For example "nearneighbor output_mean_distance.xy -R-179.5/179.5/-89.5/89.5 -I1 -N4 -S1d -Goutput_mean_distance.nc".
        gmt_command_line = [
                "gmt",
                "nearneighbor",
                xyz_filename,
                "-N4+m2", # Divide search radius into 4 sectors but only require values in 2 sectors.
                "-S{}d".format(1.5 * grid_spacing), # Search radius is a larger multiple the grid spacing.
                "-I{}".format(grid_spacing),
                # Use GMT gridline registration since our input point grid has data points on the grid lines.
                # Gridline registration is the default so we don't need to force pixel registration...
                # "-r", # Force pixel registration since data points are at centre of cells.
                "-R{}/{}/{}/{}".format(-180, 180, -90, 90),
                "-G{}".format(grd_filename)]
    else:
        # The command-line strings to execute GMT 'xyz2grd'.
        # For example "xyz2grd output_mean_distance.xy -R-179.5/179.5/-89.5/89.5 -I1 -Goutput_mean_distance.nc".
        gmt_command_line = [
                "gmt",
                "xyz2grd",
                xyz_filename,
                "-I{}".format(grid_spacing),
                # Use GMT gridline registration since our input point grid has data points on the grid lines.
                # Gridline registration is the default so we don't need to force pixel registration...
                # "-r", # Force pixel registration since data points are at centre of cells.
                "-R{}/{}/{}/{}".format(-180, 180, -90, 90),
                "-G{}".format(grd_filename)]
    
    call_system_command(gmt_command_line)

    cpu_profile.end_write_grd_file_from_xyz()


def calculate_upscaled_mask_interpolation_params(src_lon_lats, src_grid_spacing, upscaled_lon_lats_string, age_grid_filename):
    cpu_profile.start_calculate_upscaled_mask_interpolation_params()
    cpu_profile.start_upscaled_sample_age_grid()

    # Create a mask (matching age grid mask) at our upscaled grid spacing.
    #
    # Sample age grid at the upscaled grid spacing.
    # The "-s" option suppresses output of NaN values.
    upscaled_masked_lon_lat_ages_string = call_system_command(
            # The command-line strings to execute GMT 'grdtrack'...
            ["gmt", "grdtrack", "-fg", "-s", "-G{}".format(age_grid_filename)],
            stdin=upscaled_lon_lats_string,
            return_stdout=True)
    #memory_profile.print_object_memory_usage(upscaled_masked_lon_lat_ages_string, 'upscaled_masked_lon_lat_ages_string')
    
    cpu_profile.end_upscaled_sample_age_grid()
    cpu_profile.start_extract_upscale_mask_samples()

    # Extract the age-grid-masked input points.
    upscaled_masked_lon_lat_ages_lines = upscaled_masked_lon_lat_ages_string.splitlines()
    del upscaled_masked_lon_lat_ages_string  # free memory
    num_upscaled_points = len(upscaled_masked_lon_lat_ages_lines)
    upscaled_masked_lon_lats = np.empty((num_upscaled_points, 2), dtype=float)
    for line_index, line in enumerate(upscaled_masked_lon_lat_ages_lines):
        # Each line returned by GMT grdtrack contains "longitude latitude age_grid_value".
        # Note that due to "-s" option to "gmt grdtrack" we will only get non-NaN age grid values.
        lon_str, lat_str, _ = line.split()
        lon, lat = float(lon_str), float(lat_str)
        upscaled_masked_lon_lats[line_index] = (lon, lat)
    del upscaled_masked_lon_lat_ages_lines  # free memory
    #memory_profile.print_object_memory_usage(upscaled_masked_lon_lats, 'upscaled_masked_lon_lats')
    
    cpu_profile.end_extract_upscale_mask_samples()
    cpu_profile.start_create_and_query_kdtrees()

    # Create a k-d tree of the source points (in lon-lat space).
    src_kdtree = KDTree(src_lon_lats)

    cpu_profile.end_create_and_query_kdtrees()

    # The upscaled points (lon, lat) and their interpolation parameters (source indices and weights).
    # Each upscaled point has 4 indices (into the source points) and 4 associated inverse-distance weights.
    # Note: If not all 4 indices/weights are used then the unused ones get a weight of zero.
    upscaled_masked_lon_lat_index6_weight6 = np.zeros(num_upscaled_points, dtype=[('lon', 'f8'), ('lat', 'f8'), ('index', 'i4', 6), ('weight', 'f8', 6)])
    upscaled_masked_valid_index = 0

    # Query the nearest neighbours of the upscaled points in a loop.
    # This reduces memory usage quite significantly (since each loop iteration processes a subset of points and hence uses less memory).
    #
    # The number of upscaled points to query at a time (per loop iteration).
    max_upscaled_points_per_kdtree = 100*1000
    upscaled_point_base_index = 0
    while upscaled_point_base_index < num_upscaled_points:
        cpu_profile.start_create_and_query_kdtrees()

        # Create a k-d tree of the subset of upscaled points in the current loop iteration (in lon-lat space).
        upscaled_masked_kdtree = KDTree(upscaled_masked_lon_lats[upscaled_point_base_index : upscaled_point_base_index + max_upscaled_points_per_kdtree])

        # For each upscaled point search within a radius around it for source points.
        # Use a search radius that will capture at most 6 nearest neighbours (of uniformly gridded source points).
        # So we want the search radius to be between sqrt(1^2 + 1^2) = 1.414 and sqrt(1^2 + 0.5^2) = 1.118 as shown in the two diagrams...
        #
        # . x .   . x x .
        # x o x   . xox .
        # . x .   . x x .
        #
        # ...where 'x' is a near neighbour source point and 'o' is the upscaled point.
        # The first diagram has 5 source points (one 'x' is at the 'o') and the second has 6 source points.
        search_radius = 1.2 * src_grid_spacing
        upscaled_src_indices = upscaled_masked_kdtree.query_ball_tree(src_kdtree, search_radius)
        #memory_profile.print_object_memory_usage(upscaled_src_indices, 'upscaled_src_indices')

        cpu_profile.end_create_and_query_kdtrees()
        cpu_profile.start_calc_interpolation_weights()

        for upscaled_point_index, src_indices in enumerate(upscaled_src_indices):
            # It's possible there are no near neighbours within the search radius.
            # This can happen if the source points did not adequately capture long thin geographical structures in the ocean.
            if not src_indices:
                continue

            #
            # NOTE: In this loop we try NOT to call any NumPy functions/operators.
            #       The call overhead for each numpy call is around 0.25 (for multiply operator '*') to 1.5 microseconds (for np.sum).
            #       And that's before any internal calculations.
            #       Contrast that with Python's builtin 'sum' which has only 0.05 microseconds of overhead (ie, 30 times less!).
            #       At an upscale grid spacing of 0.1 degrees there are 6,485,401 points to process in this loop.
            #       So a single 'np.sum' call will add an overhead of about 10 seconds (6,485,401 * 1.5e-6).
            #       And, even worse, in the inner loop below (over nearest source points) multiply this time by about 5 (number of near neighbours)
            #       to get 50 seconds of overhead per 'np.sum' call!
            #

            upscaled_lon, upscaled_lat = upscaled_masked_lon_lats[upscaled_point_base_index + upscaled_point_index]
            # Convert np.float to native float to avoid triggering subsequent numpy calls.
            upscaled_lon, upscaled_lat = float(upscaled_lon), float(upscaled_lat)

            # Each upscaled point gets 6 source indices and 6 weights.
            upscaled_src_index6 = [0, 0, 0, 0, 0, 0]
            upscaled_src_weight6 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            # Populate up to 6 weights starting at the first element of each 6-element index/weight array.
            current_src_weight_index = 0

            # Calculate up to 6 source indices/weights surrounding the current upscaled point.
            for src_index in src_indices:
                src_lon, src_lat = src_lon_lats[src_index]
                # Convert np.float to native float to avoid triggering subsequent numpy calls.
                src_lon, src_lat =  float(src_lon), float(src_lat)

                # Distance from upscaled point to its current near neighbour.
                delta_lon = upscaled_lon - src_lon
                delta_lat = upscaled_lat - src_lat
                distance = math.sqrt(delta_lon*delta_lon + delta_lat*delta_lat)

                # Calculate the filter weight for the current source point based on its distance to the upscaled point.
                #
                # The numerator is linear from 'search_radius' at centre to ~0 at the radius.
                # The denominator is linear from 1 at centre to 3 at the radius to give the filter a steeper drop-off.
                #
                # Note: The 1.0001 multiplier ensures a slightly non-zero weight (when distance==search_radius).
                #       This ensures we never get a divide-by-zero error when normalised the weights.
                upscaled_src_weight6[current_src_weight_index] = (1.0001 * search_radius - distance) / (1 + 2 * distance / search_radius)
                upscaled_src_index6[current_src_weight_index] = src_index
                current_src_weight_index += 1

                # Make sure we don't try to use more than 6 weights.
                # Shouldn't need this if the source points are uniformly gridded at 'src_grid_spacing' and we've set the search radius correctly.
                # But check just in case (we'll be satisfied with whatever 6 weights we get).
                if current_src_weight_index == 6:
                        break
            
            # Normalise the weights.
            #
            # We should have at least one non-zero weight due to the 1.0001 multiplier above (in the weight calculation), and hence no divide-by-zero error.
            #
            # Note: The first 'current_src_weight_index' weights should be non-zero, and
            #       the remaining unused weights (if any) are zero and hence don't contribute to the sum.
            inv_sum_upscaled_src_weight6 = 1.0 / sum(upscaled_src_weight6)
            upscaled_src_weight6[0] *= inv_sum_upscaled_src_weight6
            upscaled_src_weight6[1] *= inv_sum_upscaled_src_weight6
            upscaled_src_weight6[2] *= inv_sum_upscaled_src_weight6
            upscaled_src_weight6[3] *= inv_sum_upscaled_src_weight6
            upscaled_src_weight6[4] *= inv_sum_upscaled_src_weight6
            upscaled_src_weight6[5] *= inv_sum_upscaled_src_weight6

            # Store in the output array.
            upscaled_masked_lon_lat_index6_weight6[upscaled_masked_valid_index] = (
                    upscaled_lon,
                    upscaled_lat,
                    upscaled_src_index6,
                    upscaled_src_weight6)
            upscaled_masked_valid_index += 1

        upscaled_point_base_index += max_upscaled_points_per_kdtree

        cpu_profile.end_calc_interpolation_weights()
    
    # Resize the number of upscaled points since some might not be near any source points (and hence got excluded).
    upscaled_masked_lon_lat_index6_weight6 = upscaled_masked_lon_lat_index6_weight6[:upscaled_masked_valid_index]

    # The 'copy()' avoids the above view (slice) which does not count array storage.
    #memory_profile.print_object_memory_usage(upscaled_masked_lon_lat_index6_weight6.copy(), 'upscaled_masked_lon_lat_index6_weight6')

    cpu_profile.end_calculate_upscaled_mask_interpolation_params()

    return upscaled_masked_lon_lat_index6_weight6


def write_upscaled_grd_file(grd_filename, scalars, upscaled_lon_lat_indices_weights, upscaled_grid_spacing):
    cpu_profile.start_write_upscaled_grd_file()
    cpu_profile.start_calc_upscaled_scalars()

    # Calculate the upscaled scalars using the upscaled interpolation weights.
    #
    # Note: We avoid iterating over the upscaled points (in a for loop) to avoid the numpy call overhead per loop iteration
    #       (around 0.25 to 1.5 microseconds per call) when looping over, eg, 6.5 million points (for 0.1 degree upscaled grid spacing).
    upscaled_lon_lat_scalars = np.column_stack((
            upscaled_lon_lat_indices_weights['lon'],  # lons
            upscaled_lon_lat_indices_weights['lat'],  # lats
            # For each upscaled point weight the scalars from near neighbour source points...
            np.sum(scalars[upscaled_lon_lat_indices_weights['index']] * upscaled_lon_lat_indices_weights['weight'], axis=1)))  # scalars
    #memory_profile.print_object_memory_usage(upscaled_lon_lat_scalars, 'upscaled_lon_lat_scalars')

    cpu_profile.end_calc_upscaled_scalars()

    # Convert array to a string for standard-input to GMT.
    upscaled_xyz_data = ''.join('{} {} {}\n'.format(lon, lat, scalar) for lon, lat, scalar in upscaled_lon_lat_scalars)
    #memory_profile.print_object_memory_usage(upscaled_xyz_data, 'upscaled_xyz_data')

    # The command-line strings to execute GMT 'xyz2grd'.
    call_system_command([
            "gmt",
            "xyz2grd",
            "-I{}".format(upscaled_grid_spacing),
            # Use GMT gridline registration since our input point grid has data points on the grid lines.
            # Gridline registration is the default so we don't need to force pixel registration...
            # "-r", # Force pixel registration since data points are at centre of cells.
            "-R{}/{}/{}/{}".format(-180, 180, -90, 90),
            "-G{}".format(grd_filename)],
            stdin=upscaled_xyz_data)
    
    cpu_profile.end_write_upscaled_grd_file()


# Class to hold all proximity data for a specific age grid paleo time.
class ProximityData(object):

    def __init__(self, point_lons, point_lats, output_mean_proximity, output_standard_deviation_proximity, output_proximity_with_time, clamp_mean_proximity_in_kms=None):
        self.point_lons = point_lons
        self.point_lats = point_lats

        self.output_mean_proximity = output_mean_proximity
        self.output_standard_deviation_proximity = output_standard_deviation_proximity
        self.output_proximity_with_time = output_proximity_with_time

        self.clamp_mean_proximity_in_kms = clamp_mean_proximity_in_kms

        # Accumulate proximity statistics for each ocean basin point over time.
        if self.output_mean_proximity or self.output_standard_deviation_proximity:
            # Statistics to calculate mean/standard-deviation for a single ocean basin point.
            # Each ocean basin point has (num_proximities, sum_proximities, sum_square_proximities) that start at zero.
            num_points = len(self.point_lons)
            self.num_proximities = np.zeros(num_points, dtype=float)  # numpy array uses less memory
            self.sum_proximities = np.zeros(num_points, dtype=float)  # numpy array uses less memory
            self.sum_square_proximities = np.zeros(num_points, dtype=float)  # numpy array uses less memory
            # Not all ocean basin points will necessarily have statistics (because they might have been deactivated immediately).
            self.valid_point_statistics = np.full(num_points, False, dtype=bool)  # numpy array uses less memory
        
        # Keep a record of all proximity over time.
        if self.output_proximity_with_time:
            self.time_datas = {}  # dict indexed by time
    
    def add_proximity(self, proximity_in_kms, time, ocean_basin_point_index, ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat):
        # Update the proximity statistics for the current ocean basin point.
        if self.output_mean_proximity or self.output_standard_deviation_proximity:
            self.valid_point_statistics[ocean_basin_point_index] = True
            self.num_proximities[ocean_basin_point_index] += 1
            self.sum_proximities[ocean_basin_point_index] += proximity_in_kms
            self.sum_square_proximities[ocean_basin_point_index] += proximity_in_kms * proximity_in_kms
        
        # Add proximity for the current reconstructed point to a list for the reconstruction time.
        if self.output_proximity_with_time:
            # If we haven't already, create a new array of (reconstructed_lon, reconstructed_lat, proximity_in_kms) for all points at the specified time.
            if time not in self.time_datas:
                # All points at the current 'time' are invalid unless they are assigned to.
                self.time_datas[time] = np.ma.masked_all((len(self.point_lons), 3), dtype=float)  # numpy array uses less memory
            # Store the data for the current ocean point.
            self.time_datas[time][ocean_basin_point_index] = (ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat, proximity_in_kms)
    
    # Return the list of times (added with 'add_proximity()').
    def get_times(self):
        return list(self.time_datas.keys())
    
    # Return array of proximity tuples for the specified time.
    # Actually each tuple is an array of length 3 containing (reconstructed_lon, reconstructed_lat, proximity_in_kms).
    def get_time_data(self, time):
        if self.output_proximity_with_time:
            # Return time data but remove any masked entries (where points were not added).
            return self.time_datas[time].compressed().reshape(-1, 3)  # raises KeyError if time not in dict
        else:
            return []
    
    # Return array of (lon, lat) points that have mean/standard-deviation statistics.
    # These are ocean basin points that have not been deactivated immediately.
    # Each tuple is (lon, lat).
    def get_mean_standard_deviation_lon_lats(self):
        return np.column_stack((
                self.point_lons[self.valid_point_statistics],
                self.point_lats[self.valid_point_statistics]))
    
    # Return array of means (over time).
    # Each array element is a single mean.
    # The order and number of elements is same as 'get_mean_standard_deviation_lon_lats()'.
    def get_means(self):
        if self.output_mean_proximity:
            valid_stats_mask = self.valid_point_statistics

            # Calculate a mean proximity over time for each ocean basin point (that has statistics).
            mean_proximity = self.sum_proximities[valid_stats_mask] / self.num_proximities[valid_stats_mask]

            # Clamp mean proximity if requested.
            if self.clamp_mean_proximity_in_kms is not None:
                mean_proximity[mean_proximity > self.clamp_mean_proximity_in_kms] = self.clamp_mean_proximity_in_kms
            
            return mean_proximity
        else:
            return []
    
    # Return array of standard deviations (over time).
    # Each array element is a single standard deviation.
    # The order and number of elements is same as 'get_mean_standard_deviation_lon_lats()'.
    def get_standard_deviations(self):
        if self.output_standard_deviation_proximity:
            valid_stats_mask = self.valid_point_statistics

            # Calculate a standard deviation proximity over time for each ocean basin point.
            mean_proximity = self.sum_proximities[valid_stats_mask] / self.num_proximities[valid_stats_mask]
            standard_deviation_proximity = np.sqrt(
                    (self.sum_square_proximities[valid_stats_mask] / self.num_proximities[valid_stats_mask]) - (mean_proximity * mean_proximity))

            # Ensure not negative due to numerical precision (ie, sqrt(negative_number)).
            standard_deviation_proximity[np.isnan(standard_deviation_proximity)] = 0.0
            
            return standard_deviation_proximity
        else:
            return []


def proximity(
        input_points, # List of (lon, lat) tuples.
        rotation_filenames,
        proximity_filenames,
        proximity_features_are_topological,
        proximity_feature_types,
        topological_reconstruction_filenames,
        age_grid_filenames_and_paleo_times,
        time_increment,
        output_distance_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        max_topological_reconstruction_time = None,
        continent_obstacle_filenames = None,
        anchor_plate_id = 0,
        proximity_distance_threshold_radians = None,
        clamp_mean_proximity_distance_radians = None):
    """
    Find the minimum distance of ocean basin point locations to proximity features (topological boundaries or non-topological features) over time.
    
    If an ocean basin point falls outside an age grid (in masked region) then it is ignored.

    Ocean points are not reconstructed earlier than 'max_topological_reconstruction_time'
    (each ocean point, for each age grid, is reconstructed back to its age grid value or 'max_topological_reconstruction_time', whichever is smaller).
    If it's 'None' then only the age grid limits how far back each point is reconstructed.

    Continent obstacle filenames can optionally be specified. If they are specified then they will form geometry obstacles that the
    shortest distance path must go around (ie, water must flow around continents). Obstacles can be both polygons and polylines.
    
    A threshold distance can be specified to reject proximities exceeding it.

    The age grids and their paleo times are specified with 'age_grid_filenames_and_paleo_times' which should be a sequence of 2-tuples (filename, time).
    You can specify more than one age grid, and hence generate more than one returned ProximityData object, because processing multiple age grids together
    reduces the running time compared to processing them individually. If there's not enough age grids processed together
    (and preferably grids with similar paleo times) then we'll spend too much time resolving/reconstructing proximity features
    (and generating shortest path obstacle grids). This is because each age grid involves reconstructing all its ocean points back in time until they disappear
    (at mid-ocean ridge) and so there's a lot of overlap in time across the age grids (where calculations can be shared if processed together).
    
    The proximity results are returned in as a dict mapping age grid paleo times to ProximityData objects.
    An age grid paleo time will be missing from the dict if all input points are outside the associated age grid (in masked regions).
    """
    
    # Make sure pygplates has support for TopologicalModel.
    if pygplates.Version.get_imported_version() < pygplates.Version(30):
        raise RuntimeError(
            "Using pygplates version {} but use of TopologicalModel requires version 0.30 or greater".format(
                pygplates.Version.get_imported_version()))
    
    if time_increment <= 0:
        raise ValueError('The time increment "{}" is not positive and non-zero.'.format(time_increment))
    
    age_grid_paleo_times = [age_grid_paleo_time for _, age_grid_paleo_time in age_grid_filenames_and_paleo_times]
    if any(age_grid_paleo_time < 0 for age_grid_paleo_time in age_grid_paleo_times):
        raise ValueError('Age grid paleo time must not be negative.')
    
    if (not output_distance_with_time and
        not output_mean_distance and
        not output_standard_deviation_distance):
        raise ValueError('No output specified for ocean basin proximity.')
    
    # Convert clamp-mean-proximity-distance from radians to Kms (if it was specified).
    if clamp_mean_proximity_distance_radians is not None:
        clamp_mean_proximity_distance_kms = clamp_mean_proximity_distance_radians * pygplates.Earth.mean_radius_in_kms
    else:
        clamp_mean_proximity_distance_kms = None
    
    cpu_profile.start_proximity()
    cpu_profile.start_read_input_data()
    
    rotation_model = pygplates.RotationModel(rotation_filenames, default_anchor_plate_id=anchor_plate_id)
    
    # Read/parse the proximity features once so we're not doing at each time iteration.
    proximity_features = pygplates.FeaturesFunctionArgument(proximity_filenames).get_features()

    if proximity_feature_types:
        # Create pygplates.FeatureType's from the strings.
        # We do this here since pygplates' objects are not yet pickable
        # (which is required for objects passed via multiprocessing).
        proximity_feature_types = [pygplates.FeatureType.create_from_qualified_string(feature_type)
            for feature_type in proximity_feature_types]
    
        # For *non-topological* features we can remove those not matching the allowed feature types.
        # Note that we can't do this for *topological* features because we need to resolve *all* topologies and
        # then filter out the shared topology sections by feature type.
        if not proximity_features_are_topological:
            # Create a new list containing only features matching the allowed feature types.
            proximity_features = [feature for feature in proximity_features
                    if feature.get_feature_type() in proximity_feature_types]
    
    topology_reconstruction_features = pygplates.FeaturesFunctionArgument(topological_reconstruction_filenames).get_features()
    
    if continent_obstacle_filenames:
        #print('Creating shortest path grid...')
        shortest_path_grid = shortest_path.Grid(7)  # grid spacing of ~ 0.7 degrees
        obstacle_features = pygplates.FeaturesFunctionArgument(continent_obstacle_filenames).get_features()
        #memory_profile.print_object_memory_usage(shortest_path_grid, 'shortest_path_grid')
    
    cpu_profile.end_read_input_data()
    cpu_profile.start_reconstruct_and_calculate_distances()

    # Class to manage reconstruction data for ocean basin points associated with a specific age grid / paleo time.
    class OceanBasinReconstruction(object):
        def __init__(self, lon_lat_age_list, age_grid_paleo_time):
            self.age_grid_paleo_time = age_grid_paleo_time
            self.num_points = len(lon_lat_age_list)

            # For each ocean basin point add lon-lat-point, time-of-appearance and initial-reconstructed-point to 3 separate lists.
            # The initial reconstructed point will be updated as the ocean basin points are topologically reconstructed back into time.
            # When a point is deactivated its entry is removed from 'current_reconstructed_points' and 'current_point_indices'
            # such that their lengths will decrease (possibly to zero if all points have been deactivated).
            point_lons = []
            point_lats = []
            point_ages = []
            current_point_indices = []
            current_reconstructed_points = []
            for point_index, (lon, lat, age) in enumerate(lon_lat_age_list):
                point_lons.append(lon)
                point_lats.append(lat)
                point_ages.append(age_grid_paleo_time + age)
                current_point_indices.append(point_index)
                current_reconstructed_points.append(pygplates.PointOnSphere(lat, lon))

            self.point_lons = np.array(point_lons, dtype=float)  # numpy array uses less memory
            self.point_lats = np.array(point_lats, dtype=float)  # numpy array uses less memory
            self.point_ages = np.array(point_ages, dtype=float)  # numpy array uses less memory
            self.current_point_indices = np.array(current_point_indices, dtype=int)  # numpy array uses less memory
            self.current_reconstructed_points = current_reconstructed_points
        
        def reconstruct_time_step(self, topological_model, time, time_increment):
            # Reconstruct the current points using the topological model from 'time' to 'time + time_increment'.
            # This reconstructs *backward* in time (younger to older).
            #
            # Note: Only reconstructing over a single time step at a time uses a LOT less memory than reconstructing
            #       ocean points over the full lifetime of oceanic crust (multiplied by the number of age grids).
            #       Once we extract the reconstructed points for this time step the reconstructed time span is released.
            reconstructed_time_span = topological_model.reconstruct_geometry(
                    self.current_reconstructed_points,
                    initial_time=time,
                    oldest_time=time + time_increment,
                    youngest_time=time,
                    time_increment=time_increment,
                    # Disable collision detection since currently it's creating some artefacts along topological boundaries.
                    # TODO: Improve collision detection in pyGPlates before enabling this...
                    deactivate_points=None)

            # Extract the reconstructed points at 'time + time_increment'.
            # Any deactivated points will be None.
            reconstructed_points = reconstructed_time_span.get_geometry_points(time + time_increment, return_inactive_points=True)
            
            active_reconstructed_points = []
            active_point_indices = []

            #
            # Extract reconstructed points that are still active in the topological model, and
            # remove those active points that don't exist at 'time + time_increment' according to the age grid.
            #
            if reconstructed_points:  # could be None if all points were deactivated by topological model
                for reconstructed_point_index, reconstructed_point in enumerate(reconstructed_points):
                    # Exclude reconstructed points that have been deactivated by the topological model.
                    if reconstructed_point is None:
                        continue

                    # Retire current point if the time we are reconstructing to ('time + time_increment')
                    # is older (earlier than) than the point's time of appearance (according to the age grid).
                    point_index = self.current_point_indices[reconstructed_point_index]
                    point_begin_time = self.point_ages[point_index]
                    if time + time_increment > point_begin_time:
                        continue

                    active_reconstructed_points.append(reconstructed_point)
                    active_point_indices.append(point_index)
            
            self.current_reconstructed_points = active_reconstructed_points
            self.current_point_indices = np.array(active_point_indices, dtype=int)

        def is_active(self):
            # Return True if not all points have been deactivated.
            return bool(self.current_reconstructed_points)
    
    # Dict mapping age grid paleo time to OceanBasinReconstruction.
    ocean_basin_reconstructions = {}
    
    # All proximity data to return to caller.
    # This is a dict mapping age grid paleo time to ProximityData.
    proximity_datas = {}

    # List of age grids, sorted by increasing paleo time, that we've not yet started processing/reconstructing.
    unprocessed_age_grid_filenames_and_paleo_times = sorted(age_grid_filenames_and_paleo_times, key=lambda grid_and_time: grid_and_time[1])

    # Iterate from the minimum paleo time (of all age grids) until all ocean basin point locations (for all age grids) have disappeared.
    time_index = int(math.ceil(unprocessed_age_grid_filenames_and_paleo_times[0][1] / time_increment))  # unprocessed age grids are sorted by paleo time

    while True:
        
        time = time_index * time_increment
        #print('Time {}'.format(time))

        # We cannot reconstruct further in the past than allowed by the topological reconstruction features.
        if (max_topological_reconstruction_time is not None and
            time > max_topological_reconstruction_time):
            break

        # We're creating the topological model at each time step (rather than once before all time steps) to avoid
        # the memory usage of accumulated resolved topologies (cached inside TopologicalModel) over the entire time range.
        # This makes it run a fraction slower but it's worth it to save the memory usage which is about 1GB
        # (for typical topological models), and if you have 16 CPUs running in parallel that's an extra 16GB.
        topological_model = pygplates.TopologicalModel(topology_reconstruction_features, rotation_model)

        cpu_profile.start_read_age_grid()
        
        # Add any age grids with a paleo-time younger (smaller) than the current time.
        # We'll need to start reconstructing them back through time (as ocean basin reconstructions).
        while (unprocessed_age_grid_filenames_and_paleo_times and
               time >= unprocessed_age_grid_filenames_and_paleo_times[0][1]):
            age_grid_filename, age_grid_paleo_time = unprocessed_age_grid_filenames_and_paleo_times.pop(0)

            # Get the ages of the input points.
            lon_lat_age_list = get_positions_and_ages(input_points, age_grid_filename)
            # If there are input points inside the age grid (in non-masked regions) then create an ocean basin reconstruction,
            # otherwise there will be no reconstruction associated with the current age grid.
            if lon_lat_age_list:
                ocean_basin_reconstruction = OceanBasinReconstruction(lon_lat_age_list, age_grid_paleo_time)
                del lon_lat_age_list  # free memory
                
                if ocean_basin_reconstruction.is_active():
                    # Add to the ocean basin reconstructions currently in progress.
                    ocean_basin_reconstructions[age_grid_paleo_time] = ocean_basin_reconstruction
                    # Also create a ProximityData object for the new ocean basin reconstruction.
                    proximity_datas[age_grid_paleo_time] = ProximityData(ocean_basin_reconstruction.point_lons, ocean_basin_reconstruction.point_lats,
                                                                         output_mean_distance, output_standard_deviation_distance, output_distance_with_time,
                                                                         clamp_mean_proximity_distance_kms)
                    #print('Created age grid {} at time {}'.format(age_grid_paleo_time, time))
                    memory_profile.print_object_memory_usage(ocean_basin_reconstructions[age_grid_paleo_time], 'ocean_basin_reconstructions[{}]'.format(age_grid_paleo_time))
        
        cpu_profile.end_read_age_grid()
        
        # If there are no unprocessed age grids and no associated ocean basin reconstructions in progress then we're finished.
        if (not unprocessed_age_grid_filenames_and_paleo_times and
            not ocean_basin_reconstructions):
            break
        
        cpu_profile.start_reconstruct_proximity()
        
        if proximity_features_are_topological:
            # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
            # We generate both the resolved topology boundaries and the boundary sections between them.
            proximity_resolved_topologies = []
            proximity_shared_boundary_sections = []
            pygplates.resolve_topologies(proximity_features, rotation_model, proximity_resolved_topologies, time, proximity_shared_boundary_sections)
            
            # Iterate over the shared boundary sections of all resolved topologies.
            proximity_reconstructed_geometries = []
            for proximity_shared_boundary_section in proximity_shared_boundary_sections:
                # Skip sections that are not included in the list of boundary feature types (if any).
                proximity_feature = proximity_shared_boundary_section.get_feature()
                if (proximity_feature_types and
                    proximity_feature.get_feature_type() not in proximity_feature_types):
                    continue
                
                # Iterate over the shared sub-segments of the current boundary line.
                # These are the parts of the boundary line that actually contribute to topological boundaries.
                for proximity_shared_sub_segment in proximity_shared_boundary_section.get_shared_sub_segments():
                    proximity_reconstructed_geometries.append(proximity_shared_sub_segment.get_resolved_geometry())
            
            del proximity_resolved_topologies  # free memory
            del proximity_shared_boundary_sections  # free memory
            
        else: # non-topological features...
            
            # Reconstruct the non-topological features that exist at the current 'time'.
            proximity_reconstructed_feature_geometries = []
            pygplates.reconstruct(proximity_features, rotation_model, proximity_reconstructed_feature_geometries, time)
            
            proximity_reconstructed_geometries = []
            for proximity_reconstructed_feature_geometry in proximity_reconstructed_feature_geometries:
                proximity_reconstructed_geometries.append(proximity_reconstructed_feature_geometry.get_reconstructed_geometry())

            del proximity_reconstructed_feature_geometries  # free memory
    
        cpu_profile.end_reconstruct_proximity()
        cpu_profile.start_calculate_distances()
        
        if continent_obstacle_filenames:

            cpu_profile.start_obstacle_reconstruct_resolve()

            obstacle_reconstructed_feature_geometries = []
            pygplates.reconstruct(obstacle_features, rotation_model, obstacle_reconstructed_feature_geometries, time)
            obstacle_reconstructed_geometries = [obstacle_reconstructed_feature_geometry.get_reconstructed_geometry()
                    for obstacle_reconstructed_feature_geometry in obstacle_reconstructed_feature_geometries]

            topology_obstacle_feature_types = [pygplates.FeatureType.gpml_mid_ocean_ridge, pygplates.FeatureType.gpml_subduction_zone]
            #topology_obstacle_feature_types = None
            topology_obstacle_shared_boundary_sections = topological_model.topological_snapshot(time).get_resolved_topological_sections()
            for topology_obstacle_shared_boundary_section in topology_obstacle_shared_boundary_sections:
                # Skip sections that are not included in the list of boundary feature types (if any).
                topology_obstacle_feature = topology_obstacle_shared_boundary_section.get_feature()
                if (topology_obstacle_feature_types and
                    topology_obstacle_feature.get_feature_type() not in topology_obstacle_feature_types):
                    continue
                
                for topology_obstacle_shared_sub_segment in topology_obstacle_shared_boundary_section.get_shared_sub_segments():
                    obstacle_reconstructed_geometries.append(topology_obstacle_shared_sub_segment.get_resolved_geometry())
        
            cpu_profile.end_obstacle_reconstruct_resolve()
            cpu_profile.start_create_obstacle_grids()
            
            # Create obstacle grid.
            shortest_path_obstacle_grid = shortest_path_grid.create_obstacle_grid(obstacle_reconstructed_geometries)

            # Create distance grid.
            shortest_path_distance_grid = shortest_path_obstacle_grid.create_distance_grid(
                    proximity_reconstructed_geometries, proximity_distance_threshold_radians)

            cpu_profile.end_create_obstacle_grids()
            cpu_profile.start_calculate_obstacle_distances()
            
            # Query distances to ocean points.
            # Find the shortest path distance to each the ocean basin point in each age grid currently being reconstructed (to all proximity reconstructed geometries).
            for age_grid_paleo_time, ocean_basin_reconstruction in ocean_basin_reconstructions.items():
                proximity_data = proximity_datas[age_grid_paleo_time]

                for ocean_basin_reconstructed_point_index, ocean_basin_reconstructed_point in enumerate(ocean_basin_reconstruction.current_reconstructed_points):
                    # Find minimum distance.
                    min_distance = shortest_path_distance_grid.shortest_distance(ocean_basin_reconstructed_point)
                    if min_distance is None:
                        # All proximity geometries are unreachable or further than distance threshold.
                        # Use longest great circle distance between two points on the globe to represent this.
                        min_distance = math.pi
                    distance_in_kms = min_distance * pygplates.Earth.mean_radius_in_kms

                    # Add minimum distance to proximity data.
                    ocean_basin_point_index = ocean_basin_reconstruction.current_point_indices[ocean_basin_reconstructed_point_index]
                    ocean_basin_reconstructed_lat, ocean_basin_reconstructed_lon = ocean_basin_reconstructed_point.to_lat_lon()
                    proximity_data.add_proximity(distance_in_kms, time, ocean_basin_point_index, ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat)
        
            # Remove references - might help Python to deallocate these objects now.
            #memory_profile.print_object_memory_usage(shortest_path_distance_grid, 'shortest_path_distance_grid')
            del shortest_path_distance_grid
            #memory_profile.print_object_memory_usage(shortest_path_obstacle_grid, 'shortest_path_obstacle_grid')
            del shortest_path_obstacle_grid
            del obstacle_reconstructed_feature_geometries
            del obstacle_reconstructed_geometries
            del topology_obstacle_shared_boundary_sections
            del proximity_reconstructed_geometries
        
            cpu_profile.end_calculate_obstacle_distances()
            
        else:
            # Find the minimum distance to each the ocean basin point in each age grid currently being reconstructed (to all proximity reconstructed geometries).
            for age_grid_paleo_time, ocean_basin_reconstruction in ocean_basin_reconstructions.items():

                # Find minimum distances.
                proximity_geometries_closest_to_ocean_basin_points = proximity_query.find_closest_geometries_to_points(
                        ocean_basin_reconstruction.current_reconstructed_points,
                        proximity_reconstructed_geometries,
                        distance_threshold_radians = proximity_distance_threshold_radians)
                
                # Add minimum distances to proximity data.
                proximity_data = proximity_datas[age_grid_paleo_time]
                for ocean_basin_reconstructed_point_index, proximity_geometry_closest_to_ocean_basin_point in enumerate(proximity_geometries_closest_to_ocean_basin_points):
                    if proximity_geometry_closest_to_ocean_basin_point is not None:
                        min_distance, _ = proximity_geometry_closest_to_ocean_basin_point
                    else:
                        # All proximity geometries are unreachable or further than distance threshold.
                        # Use longest great circle distance between two points on the globe to represent this.
                        min_distance = math.pi
                    distance_in_kms = min_distance * pygplates.Earth.mean_radius_in_kms

                    ocean_basin_point_index = ocean_basin_reconstruction.current_point_indices[ocean_basin_reconstructed_point_index]
                    ocean_basin_reconstructed_point = ocean_basin_reconstruction.current_reconstructed_points[ocean_basin_reconstructed_point_index]
                    ocean_basin_reconstructed_lat, ocean_basin_reconstructed_lon = ocean_basin_reconstructed_point.to_lat_lon()
                    proximity_data.add_proximity(distance_in_kms, time, ocean_basin_point_index, ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat)
                del proximity_geometries_closest_to_ocean_basin_points  # free memory
            
            del proximity_reconstructed_geometries  # free memory
    
        cpu_profile.end_calculate_distances()
        cpu_profile.start_reconstruct_time_step()
        
        # Reconstruct to the next time unless we're already at the last time.
        if (max_topological_reconstruction_time is None or
            time <= max_topological_reconstruction_time):
            
            # Find the ocean basin points for the next time step.
            # We do this for each age grid currently being reconstructed.
            for age_grid_paleo_time in list(ocean_basin_reconstructions.keys()):  # copy dict keys since might remove them while iterating
                ocean_basin_reconstruction = ocean_basin_reconstructions[age_grid_paleo_time]
                # Reconstruct the current ocean basin points from 'time' to 'time + time_increment'.
                # The reconstructed points will be the current points in the next time step.
                ocean_basin_reconstruction.reconstruct_time_step(topological_model, time, time_increment)
                # If finished reconstructing ocean basin (for associated age grid) then remove from current reconstructions.
                if not ocean_basin_reconstruction.is_active():
                    #print('Finished age grid {} at time {}'.    format(age_grid_paleo_time, time))
                    del ocean_basin_reconstructions[age_grid_paleo_time]
            
            del topological_model  # free memory
    
        cpu_profile.end_reconstruct_time_step()
        
        # Increment the time (to the next time interval).
        time_index += 1

    cpu_profile.end_reconstruct_and_calculate_distances()
    cpu_profile.end_proximity()
    
    memory_profile.print_object_memory_usage(proximity_datas, 'proximity_datas')
    #for proximity_data_time in proximity_datas.keys():
    #    memory_profile.print_object_memory_usage(proximity_datas[proximity_data_time], 'proximity_datas[{}]'.format(proximity_data_time))

    return proximity_datas
    
    
def write_proximity_data(
        proximity_datas,
        age_grid_filenames_and_paleo_times,
        output_filename_prefix,
        output_filename_extension,
        output_distance_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        output_grd_files = None):
    
    cpu_profile.start_write_proximity_data()

    # If we're outputting mean and/or standard deviation grids, and if upscaling has been enabled.
    if output_grd_files and (output_mean_distance or output_standard_deviation_distance):
        _, upscale_mean_std_dev_grid_spacing = output_grd_files
        if upscale_mean_std_dev_grid_spacing is not None:
            # Generate input points at the upscaled grid spacing.
            # These will be used to generate an upscaled mask from an age grid.
            # We do this outside the loop over age grids because it only needs to be done once (for all age grids) and so reduces running time.
            cpu_profile.start_upscaled_mask_generate_input_points()
            upscaled_lon_lats, _, _ = generate_input_points_grid(upscale_mean_std_dev_grid_spacing)  # this is quite fast (using numpy)
            upscaled_lon_lats_string = ''.join('{} {}\n'.format(lon, lat) for lon, lat in upscaled_lon_lats)  # this is quite slow
            cpu_profile.end_upscaled_mask_generate_input_points()
    
    # Write the distance grid(s) associated with each input age grid.
    for age_grid_filename, age_grid_paleo_time in age_grid_filenames_and_paleo_times:

        proximity_data = proximity_datas.get(age_grid_paleo_time)  # lookup in dict
        if not proximity_data:
            print('WARNING: All ocean basin points are outside the age grid: {}'.format(age_grid_filename), file=sys.stderr)
            continue
    
        if output_distance_with_time:

            for time in proximity_data.get_times():

                xyz_filename = '{}_{:.1f}_{:.1f}.{}'.format(output_filename_prefix, age_grid_paleo_time, time, output_filename_extension)
                write_xyz_file(xyz_filename, proximity_data.get_time_data(time))

                if output_grd_files:
                    grd_filename = '{}_{:.1f}_{:.1f}.nc'.format(output_filename_prefix, age_grid_paleo_time, time)
                    ocean_basin_grid_spacing, _ = output_grd_files
                    write_grd_file_from_xyz(
                            grd_filename, xyz_filename, ocean_basin_grid_spacing,
                            # Using reconstructed points (which are *not* grid-aligned) so need to use nearest neighbour filtering...
                            use_nearneighbor=True)

        # If we're outputting mean and/or standard deviation grids, and if upscaling has been enabled.
        if output_mean_distance or output_standard_deviation_distance:
            mean_standard_deviation_lon_lats = proximity_data.get_mean_standard_deviation_lon_lats()
            # If we're outputting mean and/or standard deviation grids, and if upscaling has been enabled.
            if output_grd_files:
                ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    upscaled_masked_lon_lat_indices_weights = calculate_upscaled_mask_interpolation_params(
                            mean_standard_deviation_lon_lats,
                            ocean_basin_grid_spacing,
                            upscaled_lon_lats_string,
                            age_grid_filename)
        
        if output_mean_distance:

            means = proximity_data.get_means()

            # Write the xyz file.
            xyz_mean_distance_filename = '{}_{:.1f}_mean_distance.{}'.format(output_filename_prefix, age_grid_paleo_time, output_filename_extension)
            # An array of (lon, lat, mean).
            xyz_mean_data = np.column_stack((mean_standard_deviation_lon_lats, means))
            write_xyz_file(xyz_mean_distance_filename, xyz_mean_data)
            
            # Write the grid file.
            if output_grd_files:
                grd_mean_distance_filename = '{}_{:.1f}_mean_distance.nc'.format(output_filename_prefix, age_grid_paleo_time)
                ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    write_upscaled_grd_file(
                            grd_mean_distance_filename,
                            means,
                            upscaled_masked_lon_lat_indices_weights,
                            upscale_mean_std_dev_grid_spacing)
                else:
                    write_grd_file_from_xyz(
                            grd_mean_distance_filename, xyz_mean_distance_filename, ocean_basin_grid_spacing,
                            # Using original (grid-aligned) points so don't near nearest neighbour filtering...
                            use_nearneighbor=False)
        
        if output_standard_deviation_distance:

            standard_deviations = proximity_data.get_standard_deviations()

            # Write the xyz file.
            xyz_standard_deviation_distance_filename = '{}_{:.1f}_std_dev_distance.{}'.format(output_filename_prefix, age_grid_paleo_time, output_filename_extension)
            # An array of (lon, lat, standard_deviation).
            xyz_standard_deviation_data = np.column_stack((mean_standard_deviation_lon_lats, standard_deviations))
            write_xyz_file(xyz_standard_deviation_distance_filename, xyz_standard_deviation_data)

            # Write the grid file.
            if output_grd_files:
                grd_standard_deviation_distance_filename = '{}_{:.1f}_std_dev_distance.nc'.format(output_filename_prefix, age_grid_paleo_time)
                ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    write_upscaled_grd_file(
                            grd_standard_deviation_distance_filename,
                            standard_deviations,
                            upscaled_masked_lon_lat_indices_weights,
                            upscale_mean_std_dev_grid_spacing)
                else:
                    write_grd_file_from_xyz(
                            grd_standard_deviation_distance_filename, xyz_standard_deviation_distance_filename, ocean_basin_grid_spacing,
                            # Using original (grid-aligned) points so don't near nearest neighbour filtering...
                            use_nearneighbor=False)
    
    cpu_profile.end_write_proximity_data()
    
    # See how much extra memory is used after time/mean/standard-deviation data is extracted from the ProximityData.
    #for proximity_data_time in proximity_datas.keys():
    #    memory_profile.print_object_memory_usage(proximity_datas[proximity_data_time], 'proximity_datas[{}] at end of write_proximity_data()'.format(proximity_data_time))


def generate_and_write_proximity_data(
        input_points, # List of (lon, lat) tuples.
        rotation_filenames,
        proximity_filenames,
        proximity_features_are_topological,
        proximity_feature_types,
        topological_reconstruction_filenames,
        age_grid_filenames_and_paleo_times,
        time_increment,
        output_distance_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        output_filename_prefix,
        output_filename_extension,
        max_topological_reconstruction_time = None,
        continent_obstacle_filenames = None,
        anchor_plate_id = 0,
        proximity_distance_threshold_radians = None,
        clamp_mean_proximity_distance_radians = None,
        output_grd_files = None):
    
    # Calculate proximity data.
    proximity_datas = proximity(
            input_points, # List of (lon, lat) tuples.
            rotation_filenames,
            proximity_filenames,
            proximity_features_are_topological,
            proximity_feature_types,
            topological_reconstruction_filenames,
            age_grid_filenames_and_paleo_times,
            time_increment,
            output_distance_with_time,
            output_mean_distance,
            output_standard_deviation_distance,
            max_topological_reconstruction_time,
            continent_obstacle_filenames,
            anchor_plate_id,
            proximity_distance_threshold_radians,
            clamp_mean_proximity_distance_radians)

    # Write proximity data.
    write_proximity_data(
            proximity_datas,
            age_grid_filenames_and_paleo_times,
            output_filename_prefix,
            output_filename_extension,
            output_distance_with_time,
            output_mean_distance,
            output_standard_deviation_distance,
            output_grd_files)
    
    # Print CPU usage.
    age_grid_paleo_times = [time for _, time in age_grid_filenames_and_paleo_times]
    cpu_profile.print_usage(age_grid_paleo_times)


# Wraps around 'generate_and_write_proximity_data()' so can be used by multiprocessing.Pool.map() which requires a single-argument function.
def generate_and_write_proximity_data_parallel_pool_function(args):
    try:
        return generate_and_write_proximity_data(*args)
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


def generate_and_write_proximity_data_parallel(
        input_points, # List of (lon, lat) tuples.
        rotation_filenames,
        proximity_filenames,
        proximity_features_are_topological,
        proximity_feature_types,
        topological_reconstruction_filenames,
        age_grid_filenames_and_paleo_times,
        time_increment,
        output_distance_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        output_filename_prefix,
        output_filename_extension,
        max_topological_reconstruction_time = None,
        continent_obstacle_filenames = None,
        anchor_plate_id = 0,
        proximity_distance_threshold_radians = None,
        clamp_mean_proximity_distance_radians = None,
        output_grd_files = None,
        num_cpus = None,  # if None then defaults to all available CPUs
        max_memory_usage_in_gb = None):  # max memory to use (in GB)
    
    # If the user requested all available CPUs then attempt to find out how many there are.
    if not num_cpus:
        try:
            num_cpus = multiprocessing.cpu_count()
        except NotImplementedError:
            num_cpus = 1
    
    num_age_grids = len(age_grid_filenames_and_paleo_times)

    # Give each task a reasonable number of age grids (times) to process - if there's not enough times per task then we'll
    # spend too much time resolving/reconstructing proximity features (and generating shortest path obstacle grids) -
    # which needs to be repeated for each task (group of times) - this is because each age grid involves
    # reconstructing all its ocean points back in time until they disappear (at mid-ocean ridge) and so there's
    # a lot of overlap in time across the age grids (where calculations can be shared within a single process).
    min_num_age_grids_per_task = 10
    
    # These are rough figures determined empirically by running this script (and using MemoryProfile).
    #
    # The base amount of memory usage per task (in GB) to set up for processing (excludes age-grid/ocean-basin reconstructions).
    #
    # Note: Each increment/decrement of subdivision depth in 'shortest_path.Grid' increases/decreases usage by ~0.5 GB.
    base_memory_usage_per_task_in_gb = 2.8
    # The memory usage per age grid is roughly proportional to the number of input points,
    # with a uniform lon-lat grid at 1 degree resolution consuming about 6MB.
    delta_memory_usage_per_age_grid_in_gb = 6e-3 * len(input_points) / (180 * 360)
    # The total memory used to process the specified number of age grids in a single task.
    def memory_usage_per_task(num_age_grids_per_task_):
        return base_memory_usage_per_task_in_gb + num_age_grids_per_task_ * delta_memory_usage_per_age_grid_in_gb

    # Did we clamp to the minimum number of age grids per task?
    clamped_minimum_num_age_grids_per_task = False

    # If we've been given a limit on memory usage then determine how many age grids to process per task.
    if max_memory_usage_in_gb:
        # The memory used by the number of age grids per task multiplied by the number of tasks processed in parallel should not exceed the maximum memory usage.
        num_age_grids_per_task = math.trunc(((max_memory_usage_in_gb / num_cpus) - base_memory_usage_per_task_in_gb) / delta_memory_usage_per_age_grid_in_gb)  # could be negative
        # But don't reduce below the minimum number of age grids per task.
        if num_age_grids_per_task < min_num_age_grids_per_task:
            num_age_grids_per_task = min_num_age_grids_per_task
            clamped_minimum_num_age_grids_per_task = True
            # Reduce the number of CPUs to compensate for the higher than expected number of age grids per task, so that we don't exceed the memory limit.
            # Number of CPUs is the max memory divided by the memory used to process 'num_age_grids_per_task' age grids.
            num_cpus = math.trunc(max_memory_usage_in_gb / memory_usage_per_task(num_age_grids_per_task))
            if num_cpus < 1:
                num_cpus = 1
    else:
        # No limits on memory usage were specified, so just use the minimum number of age grids per task.
        num_age_grids_per_task = min_num_age_grids_per_task
    
    num_total_tasks = math.ceil(num_age_grids / num_age_grids_per_task)

    # If there are fewer total tasks than twice the number of CPUs then reduce the number of age grids per task
    # (but not less than the minimum) so that each CPU gets two tasks to process.
    # This helps to better utilise all CPUs when each task takes a different amount of time to complete.
    # But this is only useful if using more than one CPU.
    if num_cpus > 1:
        if num_total_tasks < 2 * num_cpus:
            # Note that reducing the number of age grids per task uses less memory (ie, doesn't violate our max memory limit).
            num_age_grids_per_task = math.ceil(num_age_grids / (2 * num_cpus))
            clamped_minimum_num_age_grids_per_task = False  # just re-calculated 'num_age_grids_per_task'
            if num_age_grids_per_task < min_num_age_grids_per_task:
                num_age_grids_per_task = min_num_age_grids_per_task
                clamped_minimum_num_age_grids_per_task = True
            num_total_tasks = math.ceil(num_age_grids / num_age_grids_per_task)
            # If there are fewer tasks than the number of CPUs reduce the number of CPUs to match.
            if num_total_tasks < num_cpus:
                num_cpus = num_total_tasks
    
    # Number of age grids per task does not need to exceed the total number of age grids.
    # This happens when there's just one task.
    if num_age_grids_per_task > num_age_grids:
        num_age_grids_per_task = num_age_grids

    if max_memory_usage_in_gb:
        print('Maximum memory usage: {:.2f}GB'.format(max_memory_usage_in_gb))
    print('Approximate memory usage: {:.2f}GB'.format(num_cpus * memory_usage_per_task(num_age_grids_per_task)))
    print('Number of age grids: {}'.format(num_age_grids))
    print('Number of age grids per task: {}{}'.format(num_age_grids_per_task, ' (clamped to minimum)' if clamped_minimum_num_age_grids_per_task else ''))
    print('Number of total tasks: {}'.format(num_total_tasks))
    print('Number of tasks in a batch (num CPUs): {}'.format(num_cpus))
    print('Number of task batches: {}'.format(math.ceil(num_total_tasks / num_cpus)))

    # Create a list of time periods.
    # Each task will be passed the times within a single time period.
    task_age_grid_filenames_and_paleo_times_lists = []
    task_start_time_index = 0
    while task_start_time_index < num_age_grids:
        # Pass 'num_age_grids_per_task' consecutive times into each task since the distance calculations are more efficient that way.
        task_age_grid_filenames_and_paleo_times_lists.append(
                age_grid_filenames_and_paleo_times[task_start_time_index : task_start_time_index + num_age_grids_per_task])
        task_start_time_index += num_age_grids_per_task
    
    #
    # No need for parallelisation if number of CPUs is one.
    #
    # Also can use this when there are exceptions in order to determine which source code line.
    # Because once goes through multiprocessing pools then lose error locations in source code.
    #
    if num_cpus == 1:
        for task_age_grid_filenames_and_paleo_times_list in task_age_grid_filenames_and_paleo_times_lists:
            generate_and_write_proximity_data(
                    input_points, # List of (lon, lat) tuples.
                    rotation_filenames,
                    proximity_filenames,
                    proximity_features_are_topological,
                    proximity_feature_types,
                    topological_reconstruction_filenames,
                    task_age_grid_filenames_and_paleo_times_list,
                    time_increment,
                    output_distance_with_time,
                    output_mean_distance,
                    output_standard_deviation_distance,
                    output_filename_prefix,
                    output_filename_extension,
                    max_topological_reconstruction_time,
                    continent_obstacle_filenames,
                    anchor_plate_id,
                    proximity_distance_threshold_radians,
                    clamp_mean_proximity_distance_radians,
                    output_grd_files)
        return
    
    # Split the workload across the CPUs.
    try:
        pool = multiprocessing.Pool(num_cpus, initializer=low_priority)
        pool_map_async_result = pool.map_async(
                generate_and_write_proximity_data_parallel_pool_function,
                (
                    (
                        input_points,
                        rotation_filenames,
                        proximity_filenames,
                        proximity_features_are_topological,
                        proximity_feature_types,
                        topological_reconstruction_filenames,
                        task_age_grid_filenames_and_paleo_times_list,
                        time_increment,
                        output_distance_with_time,
                        output_mean_distance,
                        output_standard_deviation_distance,
                        output_filename_prefix,
                        output_filename_extension,
                        max_topological_reconstruction_time,
                        continent_obstacle_filenames,
                        anchor_plate_id,
                        proximity_distance_threshold_radians,
                        clamp_mean_proximity_distance_radians,
                        output_grd_files
                    ) for task_age_grid_filenames_and_paleo_times_list in task_age_grid_filenames_and_paleo_times_lists
                ),
                1) # chunksize
        
        # Apparently if we use pool.map_async instead of pool.map and then get the results
        # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
        # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        try:
            pool_map_async_result.get(999999)
        except KeyboardInterrupt:
            # Note: 'finally' block below gets executed before returning.
            return
    finally:
        pool.close()
        pool.join()


if __name__ == '__main__':
    
    import traceback
    
    __description__ = \
    """Find the minimum distance of ocean basin point locations to proximity features (topological boundaries or non-topological features) over time.
    
    An input xy file can be specified. If one is not specified then a uniform lon/lat grid of points will be generated internally
    (its grid spacing is controlled with the '-i' option and the grid point latitudes/longitudes are offset by half the grid spacing,
    eg, for a 1 degree spacing the latitudes are -89.5, -88.5, ..., 89.5). If an input xy file is specified then it can contain
    arbitrary (lon, lat) point locations (it can also contain extra 3rd, 4th, etc, columns but they are ignored).
    
    The proximity features are assumed topological by default (specify '-n' if non-topological). All features
    found for non-topological files, or all resolved boundary section features for topological boundary files,
    are tested for proximity to ocean basin points by default. To restrict to specific feature types
    specify the '-b' option (eg, for proximity to subduction zones specify '-b SubductionZone').

    Optional continent obstacles can be specified that the shortest distance path must go around (ie, water flowing around continents, rather than through).
    If specified then mid-ocean ridges and subduction zones of resolved topologies are also added as obstacles (that water, and hence sediment, cannot pass through).
    If *not* specifed then distances are minimum straight-line (great circle arc) distances from ocean points to proximity geometries.
    
    The output (for each input age grid) can be any or all of the following:
     1) For each time from the paleo time of the age grid to the maximum age of all ocean basin points (determined by the age grid) an output xyz file is
        generated containing the reconstructed (lon, lat) point locations (as x,y) and the minimum distance to all proximity features (as z).
     2) A single output xyz file is generated containing (lon, lat) ocean basin point locations (as x,y) (at the paleo time of the age grid) and
        the mean distances of each point to all proximity features averaged over the point's lifetime (as determined by the age grid).
        See the '-j' option. A similar output xyz file can be generated containing standard deviation (instead of mean) distance. See the '-k' option.
    And a GMT grd file can be generated for each xyz file. See the '-w' option.
    
    If an ocean basin point falls outside an age grid then it is ignored.
    
    A threshold distance can be specified to reject proximities exceeding it. See the '-q' option.

    Note that you can specify more than one age grid and hence generate more than one mean-distance grid, for example.
    This is because processing multiple age grids together reduces the running time compared to processing them individually.
    If there's not enough age grids processed together (and preferably grids with similar paleo times) then we'll
    spend too much time resolving/reconstructing proximity features (and generating shortest path obstacle grids).
    This is because each age grid involves reconstructing all its ocean points back in time until they disappear (at mid-ocean ridge)
    and so there's a lot of overlap in time across the age grids (where calculations can be shared if processed together).

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -m topologies.gpml -b SubductionZone -s static_polygons.gpml -g age_grid_0.nc 0 age_grid_1.nc 1 -d -w -j -k -- input.xy output
     """
    
    try:
        # The command-line parser.
        parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
                metavar='rotation_filename', help='One or more rotation files.')
        parser.add_argument('-a', '--anchor', type=int, default=0,
                dest='anchor_plate_id',
                help='Anchor plate id used for reconstructing. Defaults to zero.')
        parser.add_argument('-m', '--proximity_filenames', type=str, nargs='+', required=True,
                metavar='proximity_filename',
                help='One or more proximity files to test proximity to the ocean basin points. '
                     'These can be topological or non-topological.')
        parser.add_argument('-n', '--non_topological_proximity_features', action='store_true',
                help='Whether the proximity features specified with (-m/--proximity_filenames) are non-topological. '
                     'By default they are topological.')
        parser.add_argument('-b', '--proximity_feature_types', type=str, nargs='+',
                metavar='proximity_feature_type',
                help='The feature type(s) to select proximity features. '
                     'The format should match the format of '
                     'http://www.gplates.org/docs/pygplates/generated/pygplates.FeatureType.html#pygplates.FeatureType.get_name . '
                     'For example, proximity to subduction zones is specified as SubductionZone (without the gpml: prefix). '
                     'Defaults to all proximity features (all resolved subsegments for topological features '
                     'or all proximity reconstructed feature geometries for non-topological features).')
        parser.add_argument('-o', '--continent_obstacle_filenames', type=str, nargs='+',
                metavar='continent_obstacle_filename',
                help='Optional continent obstacles that the shortest distance path must go around (ie, water flowing around continents, rather than through). '
                     'If specified then mid-ocean ridges and subduction zones of resolved topologies are also added as obstacles (that water, and hence sediment, cannot pass through). '
                     'If not specifed then distances are minimum straight-line (great circle arc) distances from ocean points to proximity geometries. '
                     'Obstacles can be both polygons and polylines. Default is no obstacles.')
        parser.add_argument('-s', '--topological_reconstruction_filenames', type=str, nargs='+', required=True,
                metavar='topological_reconstruction_filename',
                help='The filenames of the topological files used to incrementally reconstruct the (paleo) ocean basin points.')
        parser.add_argument('-x', '--max_topological_reconstruction_time', type=int,
                metavar='max_topological_reconstruction_time',
                help='Optional maximum (largest) time that the topological reconstruction files can reconstruct back to (in Ma). '
                     'Ocean points are not reconstructed earlier than this time '
                     '(each ocean point is reconstructed back to its age grid value or this value, whichever is smaller). '
                     'Value must be an integer. If not specified then only the age grid limits how far back each point is reconstructed.')
        parser.add_argument('-g', '--age_grid_filenames_format', type=str, required=True,
                help='The format string to generate age grid filenames (using the age grid paleo times). '
                     'For example, "/Users/me/AgeGridData/Muller2019-Young2019-Cao2020_AgeGrid-{:.0f}.nc" where "{:.0f}" (see Python\'s str.format() function) '
                     'will get replaced with each age grid time (from the "--age_grid_paleo_times" option) using zero decimal places (ie, integer times).')
        parser.add_argument('-p', '--age_grid_paleo_times', type=float, nargs='+', required=True,
                help='The age grid paleo times. These will be used together with the age grid filenames format (from the "--age_grid_filenames_format" option) to generate '
                     'the age grid filenames. Note that it is more efficient to process multiple age grids together, and best if their times are close together.')
        parser.add_argument('-t', '--time_increment', type=int, default=1,
                help='The time increment in My. Value must be an integer. Defaults to 1 My.')
        parser.add_argument('-q', '--max_distance_threshold', type=float,
                help='Distances above this optional maximum distance threshold (in Kms) are *ignored*. '
                     'If specified then distances (between ocean basin points and proximity features) exceeding this threshold will be *ignored* '
                     '(ocean basin point will not get output or will not contribute to mean / standard deviation), otherwise all distances are included.')
        parser.add_argument('--clamp_mean_distance', type=float,
                help='*Mean* distances (in Kms) above this optional maximum mean distance are *clamped* to it. '
                     'If specified then *mean* distances (between ocean basin points and proximity features) exceeding this value will be *clamped* to it, '
                     'otherwise mean distances are unclamped.')
        parser.add_argument('-c', '--num_cpus', type=int,
                help='The number of CPUs to use for calculations. Defaults to all available CPUs.')
        parser.add_argument('-u', '--max_memory_usage', type=int,
                dest='max_memory_usage_in_gb',
                help='The maximum amount of memory (in GB) to use (divided across the CPUs). '
                     'Should ideally be set to the amount of physical RAM (or less). Defaults to unlimited.')
        
        parser.add_argument('-d', '--output_distance_with_time', action='store_true',
                help='For each input point at each time during its lifetime write its distance to the nearest feature. '
                     'If no output options are specified then this one is used.')
        parser.add_argument('-j', '--output_mean_distance', action='store_true',
                help='For each input point write its mean distance to features averaged over its lifetime. '
                     'By default it is not written.')
        parser.add_argument('-k', '--output_std_dev_distance', action='store_true',
                help='For each input point write its standard deviation of distances to features averaged over its lifetime. '
                     'By default it is not written.')
        parser.add_argument('-w', '--output_grd_files', action='store_true',
                help='Also generate a grd file (".nc") for each xyz file. '
                     'By default only xyz files are written. '
                     'Can only be specified if "ocean_basin_points_filename" is not specified '
                     '(ie, ocean basin points must be on a uniform lon/lat grid).')
        
        parser.add_argument('-i', '--ocean_basin_grid_spacing', type=float,
                help='The grid spacing (in degrees) of ocean basin points in lon/lat space. '
                     'The grid points follow GMT gridline registration (ie, include the poles, and start/end on the dateline). '
                     'Can only be specified if "ocean_basin_points_filename" is not specified.')
        parser.add_argument('--upscale_mean_std_dev_grid_spacing', type=float,
                help='The grid spacing (in degrees) of mean and standard deviation distance grids can optionally be upscaled '
                     '(from the grid spacing of "--ocean_basin_grid_spacing"). '
                     'The grid points follow GMT gridline registration (ie, include the poles, and start/end on the dateline). '
                     'Can only be specified if "output_grd_files" is also specified but not "ocean_basin_points_filename".')
        
        def parse_unicode(value_string):
            return value_string
        
        parser.add_argument('ocean_basin_points_filename', type=parse_unicode, nargs='?',
                metavar='ocean_basin_points_filename',
                help='Optional input xy file containing the ocean basin point locations. '
                     'If not specified then a uniform lon/lat grid of points is generated. '
                     'Can only be specified if "ocean_basin_grid_spacing", "upscale_mean_std_dev_grid_spacing" '
                     'and "output_grd_files" are not specified.')
        
        parser.add_argument('output_filename_prefix', type=parse_unicode,
                metavar='output_filename_prefix',
                help='The output xy filename prefix used in all output filenames.')
        parser.add_argument('-e', '--output_filename_extension', type=str, default='xy',
                metavar='output_filename_extension',
                help='The output xy filename extension. Defaults to "xy".')
        
        # Parse command-line options.
        args = parser.parse_args()
        
        # Default to 'output_distance_with_time' if no output options are specified.
        if (not args.output_distance_with_time and
            not args.output_mean_distance and
            not args.output_std_dev_distance):
            args.output_distance_with_time = True
        
        # Generate a list of tuples of age grid filename and paleo time.
        # These are generated from the age grid filenames format and a list of age grid paleo times.
        age_grid_filenames_and_paleo_times = [(args.age_grid_filenames_format.format(age_grid_paleo_time), age_grid_paleo_time)
                                              for age_grid_paleo_time in args.age_grid_paleo_times]
        
        if args.ocean_basin_points_filename is not None:
            if args.ocean_basin_grid_spacing is not None:
                raise argparse.ArgumentTypeError("'ocean_basin_grid_spacing' and 'ocean_basin_points_filename' cannot both be specified.")
            if args.upscale_mean_std_dev_grid_spacing is not None:
                raise argparse.ArgumentTypeError("'upscale_mean_std_dev_grid_spacing' and 'ocean_basin_points_filename' cannot both be specified.")
            if args.output_grd_files is not None:
                raise argparse.ArgumentTypeError("'output_grd_files' and 'ocean_basin_points_filename' cannot both be specified.")
        else:  # args.ocean_basin_points_filename not specified...
            if args.ocean_basin_grid_spacing is None:
                raise argparse.ArgumentTypeError("'ocean_basin_grid_spacing' must be specified if 'ocean_basin_points_filename' is not specified.")
            if args.upscale_mean_std_dev_grid_spacing is not None:
                if args.output_grd_files is None:
                    raise argparse.ArgumentTypeError("'upscale_mean_std_dev_grid_spacing' can only be specified if 'output_grd_files' is also specified.")
        
        # Get the input points.
        if args.ocean_basin_points_filename is not None:
            input_points = read_input_points(args.ocean_basin_points_filename)
        else:
            input_points, num_grid_longitudes, num_grid_latitudes = generate_input_points_grid(args.ocean_basin_grid_spacing)
        
        # Convert maximum proximity distance from Kms to radians (if it was specified).
        if args.max_distance_threshold is None:
            proximity_distance_threshold_radians = None
        else:
            proximity_distance_threshold_radians = args.max_distance_threshold / pygplates.Earth.mean_radius_in_kms
            if proximity_distance_threshold_radians > 2 * math.pi:
                # Exceeds circumference of Earth so no need for threshold.
                proximity_distance_threshold_radians = None
        
        # Convert clamp mean proximity distance from Kms to radians (if it was specified).
        if args.clamp_mean_distance is None:
            clamp_mean_proximity_distance_radians = None
        else:
            clamp_mean_proximity_distance_radians = args.clamp_mean_distance / pygplates.Earth.mean_radius_in_kms
            if clamp_mean_proximity_distance_radians > 2 * math.pi:
                # Exceeds circumference of Earth so no need for clamping.
                clamp_mean_proximity_distance_radians = None
        
        generate_and_write_proximity_data_parallel(
                input_points,
                args.rotation_filenames,
                args.proximity_filenames,
                not args.non_topological_proximity_features, # proximity_features_are_topological
                args.proximity_feature_types,
                args.topological_reconstruction_filenames,
                age_grid_filenames_and_paleo_times,
                args.time_increment,
                args.output_distance_with_time,
                args.output_mean_distance,
                args.output_std_dev_distance,
                args.output_filename_prefix,
                args.output_filename_extension,
                args.max_topological_reconstruction_time,
                args.continent_obstacle_filenames,
                args.anchor_plate_id,
                proximity_distance_threshold_radians,
                clamp_mean_proximity_distance_radians,
                (args.ocean_basin_grid_spacing, args.upscale_mean_std_dev_grid_spacing) if args.output_grd_files else None,
                args.num_cpus,
                args.max_memory_usage_in_gb)
        
        sys.exit(0)
    
    except KeyboardInterrupt:
        sys.exit(1)
    except Exception as exc:
        print('ERROR: {}'.format(exc), file=sys.stderr)
        # Uncomment this to print traceback to location of raised exception.
        #traceback.print_exc()
        
        sys.exit(1)
