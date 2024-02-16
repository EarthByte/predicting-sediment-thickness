
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
    import ptt.utils.points_in_polygons as points_in_polygons
    import ptt.utils.proximity_query as proximity_query
except ImportError:
    from gplately.ptt.utils.call_system_command import call_system_command
    import gplately.ptt.utils.points_in_polygons as points_in_polygons
    import gplately.ptt.utils.proximity_query as proximity_query
import pygplates
import shortest_path
import sys
import tempfile
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
            self.time_usage_topology_resolve_time_step = 0.0
            self.time_usage_topology_reconstruct_time_step = 0.0
            self.time_usage_topology_reconstruct_time_step_deactivate = 0.0
            self.time_usage_topology_reconstruct_time_step_find_polygons = 0.0
            self.time_usage_topology_reconstruct_time_step_stage_rotations = 0.0
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
    
    def start_topology_resolve_time_step(self):
        if self.enable_profiling:
            self.time_snapshot_start_topology_resolve_time_step = self.profile()
    def end_topology_resolve_time_step(self):
        if self.enable_profiling:
            self.time_usage_topology_resolve_time_step += self.profile() - self.time_snapshot_start_topology_resolve_time_step
    
    def start_topology_reconstruct_time_step(self):
        """Call at the start of topology_reconstruct_time_step()."""
        if self.enable_profiling:
            self.time_snapshot_start_topology_reconstruct_time_step = self.profile()
    def end_topology_reconstruct_time_step(self):
        """Call at the end of topology_reconstruct_time_step()."""
        if self.enable_profiling:
            self.time_usage_topology_reconstruct_time_step += self.profile() - self.time_snapshot_start_topology_reconstruct_time_step
    
    def start_topology_reconstruct_time_step_deactivate(self):
        if self.enable_profiling:
            self.time_snapshot_start_topology_reconstruct_time_step_deactivate = self.profile()
    def end_topology_reconstruct_time_step_deactivate(self):
        if self.enable_profiling:
            self.time_usage_topology_reconstruct_time_step_deactivate += self.profile() - self.time_snapshot_start_topology_reconstruct_time_step_deactivate
    
    def start_topology_reconstruct_time_step_find_polygons(self):
        if self.enable_profiling:
            self.time_snapshot_start_topology_reconstruct_time_step_find_polygons = self.profile()
    def end_topology_reconstruct_time_step_find_polygons(self):
        if self.enable_profiling:
            self.time_usage_topology_reconstruct_time_step_find_polygons += self.profile() - self.time_snapshot_start_topology_reconstruct_time_step_find_polygons
    
    def start_topology_reconstruct_time_step_stage_rotations(self):
        if self.enable_profiling:
            self.time_snapshot_start_topology_reconstruct_time_step_stage_rotations = self.profile()
    def end_topology_reconstruct_time_step_stage_rotations(self):
        if self.enable_profiling:
            self.time_usage_topology_reconstruct_time_step_stage_rotations += self.profile() - self.time_snapshot_start_topology_reconstruct_time_step_stage_rotations
    
    def print_proximity_usage(self, age_grid_paleo_times):
        """Call at the end of proximity() to print its CPU usage."""
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
            print(f"          Resolve topologies: {self.time_usage_topology_resolve_time_step * scale_to_seconds:.2f} seconds")
            print(f"          Topology reconstruct time step: {self.time_usage_topology_reconstruct_time_step * scale_to_seconds:.2f} seconds")
            print(f"            Topology reconstruct time step (deactivate): {self.time_usage_topology_reconstruct_time_step_deactivate * scale_to_seconds:.2f} seconds")
            print(f"            Topology reconstruct time step (points in polygons): {self.time_usage_topology_reconstruct_time_step_find_polygons * scale_to_seconds:.2f} seconds")
            print(f"            Topology reconstruct time step (stage rotations): {self.time_usage_topology_reconstruct_time_step_stage_rotations * scale_to_seconds:.2f} seconds")
    
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

    input_points = np.full((num_latitudes * num_longitudes, 2), (0.0, 0.0), dtype=float)  # numpy array uses less memory
    input_point_index = 0
    for lat_index in range(num_latitudes):
        lat = -90 + lat_index * grid_spacing_degrees
        
        for lon_index in range(num_longitudes):
            lon = -180 + lon_index * grid_spacing_degrees
            
            input_points[input_point_index] = (lon, lat)
            input_point_index += 1

    # num_latitudes = int(math.floor(180.0 / grid_spacing_degrees))
    # num_longitudes = int(math.floor(360.0 / grid_spacing_degrees))
    # for lat_index in range(num_latitudes):
    #     # The 0.5 puts the point in the centre of the grid pixel.
    #     # This also avoids sampling right on the poles.
    #     lat = -90 + (lat_index + 0.5) * grid_spacing_degrees
    #     
    #     for lon_index in range(num_longitudes):
    #         # The 0.5 puts the point in the centre of the grid pixel.
    #         # This also avoids sampling right on the dateline where there might be
    #         # age grid or static polygon artifacts.
    #         lon = -180 + (lon_index + 0.5) * grid_spacing_degrees
    #         
    #         input_points.append((lon, lat))
    
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


# Reconstruct points from 'time' to 'time + time_increment' using topologies.
def topology_reconstruct_time_step(
        ocean_basin_reconstruction,
        time,
        time_increment,
        resolved_plate_boundaries,
        resolved_plate_polygons,
        rotation_model):
    #
    # The following code comment is the *slower* version of point-in-polygon testing.
    #
    
    #    next_point_infos = []
    #
    #    for point_feature, point_begin_time, curr_recon_point in curr_point_infos:
    #        plate_id = None
    #        for resolved_plate_boundary in resolved_plate_boundaries:
    #            resolved_plate_polygon = resolved_plate_boundary.get_resolved_boundary()
    #            if resolved_plate_polygon.is_point_in_polygon(curr_recon_point):
    #                plate_id = resolved_plate_boundary.get_feature().get_reconstruction_plate_id()
    #                break
    #         
    #        if plate_id is None:
    #            # Retire current point if falls outside a topological boundary, since we don't know how to move it.
    #            # Shouldn't really happen unless falls in a gap somewhere in global topological plates.
    #            continue
    #        
    #        # Get the stage rotation from 'time' to 'time + time_increment'.
    #        stage_rotation = rotation_model.get_rotation(time + time_increment, plate_id, time)
    #        next_recon_point = stage_rotation * curr_recon_point
    #        
    #        next_point_infos.append((point_feature, point_begin_time, next_recon_point))
    #
    #    return next_point_infos
    
    #
    # The following code is the *faster* version of the above point-in-polygon code.
    #

    cpu_profile.start_topology_reconstruct_time_step()
    cpu_profile.start_topology_reconstruct_time_step_deactivate()
    
    # Read OceanBasinReconstruction member variables into local variables.
    point_ages = ocean_basin_reconstruction.point_ages
    current_reconstructed_points = ocean_basin_reconstruction.current_reconstructed_points
    current_point_indices = ocean_basin_reconstruction.current_point_indices

    # Remove any reconstructed points that don't exist at 'time + time_increment'.
    next_reconstructed_points = []
    next_point_indices = []
    for reconstructed_point_index, point_index in enumerate(current_point_indices):
        # Retire current point if the time we are reconstructing to ('time + time_increment')
        # is older (earlier than) than the point's time of appearance.
        point_begin_time = point_ages[point_index]
        if time + time_increment > point_begin_time:
            continue

        next_reconstructed_points.append(current_reconstructed_points[reconstructed_point_index])
        next_point_indices.append(point_index)

    current_reconstructed_points = next_reconstructed_points
    current_point_indices = np.array(next_point_indices, dtype=int)  # numpy array uses less memory

    cpu_profile.end_topology_reconstruct_time_step_deactivate()
    
    # If their are any active points then reconstruct them to the next time step.
    if current_reconstructed_points:

        cpu_profile.start_topology_reconstruct_time_step_find_polygons()
        
        # Find the polygon plates containing the points.
        resolved_plate_boundaries_containing_points = points_in_polygons.find_polygons(
                current_reconstructed_points,
                resolved_plate_polygons,
                resolved_plate_boundaries)

        cpu_profile.end_topology_reconstruct_time_step_find_polygons()
        cpu_profile.start_topology_reconstruct_time_step_stage_rotations()

        # Cache the stage rotation associated with each resolved plate boundary.
        stage_rotation_cache = {}

        # Iterate over the point-in-polygon results.
        for reconstructed_point_index, resolved_plate_boundary in enumerate(resolved_plate_boundaries_containing_points):
            # If point is not within any resolved plate boundaries then don't move it.
            if resolved_plate_boundary is None:
                continue
            
            # Get the stage rotation from 'time' to 'time + time_increment'.
            if resolved_plate_boundary in stage_rotation_cache:
                stage_rotation = stage_rotation_cache[resolved_plate_boundary]
            else:
                # Only call these functions for a new resolved plate boundary (to reduce computations).
                plate_id = resolved_plate_boundary.get_feature().get_reconstruction_plate_id()
                stage_rotation = rotation_model.get_rotation(time + time_increment, plate_id, time)
                # Cache stage rotation.
                stage_rotation_cache[resolved_plate_boundary] = stage_rotation
            
            # Update the current reconstructed point position.
            current_reconstructed_points[reconstructed_point_index] = stage_rotation * current_reconstructed_points[reconstructed_point_index]

        cpu_profile.end_topology_reconstruct_time_step_stage_rotations()

    # Write modified OceanBasinReconstruction member variables from local variables.
    ocean_basin_reconstruction.current_reconstructed_points = current_reconstructed_points
    ocean_basin_reconstruction.current_point_indices = current_point_indices
    
    cpu_profile.end_topology_reconstruct_time_step()


def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def write_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, use_nearneighbor = True):
    
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


def get_upscaled_mask_grd_file(upscaled_grid_spacing, age_grid_filename):
    
    #tprof_start = time_profile.perf_counter()

    # Temporary upscaled mask grid file.
    upscaled_mask_grd_file = tempfile.NamedTemporaryFile(delete=True)
    upscaled_mask_grd_file.close()  # cannot open twice on Windows - close before opening again

    #tprof_upscaled_input_points_start = time_profile.perf_counter()

    # Generate input points at the upscaled grid spacing.
    upscaled_input_points, _, _ = generate_input_points_grid(upscaled_grid_spacing)
    upscaled_input_points_data = ''.join('{} {}\n'.format(lon, lat) for lon, lat in upscaled_input_points)

    #tprof_upscaled_input_points_end = time_profile.perf_counter()
    #print(f"  generate upscaled input points: {tprof_upscaled_input_points_end - tprof_upscaled_input_points_start:.2f} seconds")
    #tprof_upscaled_mask_data_start = time_profile.perf_counter()

    # Sample age grid at the upscaled grid spacing.
    upscaled_mask_data = call_system_command(
            # The command-line strings to execute GMT 'grdtrack'...
            ["gmt", "grdtrack", "-fg", "-G{}".format(age_grid_filename)],
            stdin=upscaled_input_points_data,
            return_stdout=True)

    # Write age grid (at upscaled grid spacing) to the upscaled mask grid file.
    call_system_command([
            "gmt",
            "xyz2grd",
            "-I{}".format(upscaled_grid_spacing),
            # Use GMT gridline registration since our input point grid has data points on the grid lines.
            # Gridline registration is the default so we don't need to force pixel registration...
            # "-r", # Force pixel registration since data points are at centre of cells.
            "-R{}/{}/{}/{}".format(-180, 180, -90, 90),
            "-G{}".format(upscaled_mask_grd_file.name)],
            stdin=upscaled_mask_data)
    
    #tprof_upscaled_mask_data_end = time_profile.perf_counter()
    #print(f"  create mask: {tprof_upscaled_mask_data_end - tprof_upscaled_mask_data_start:.2f} seconds")
    #tprof_end = time_profile.perf_counter()
    #print(f"get_upscaled_mask_grd_file: {tprof_end - tprof_start:.2f} seconds")

    return upscaled_mask_grd_file


def write_upscaled_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, upscaled_grid_spacing, upscaled_mask_grd_filename):
    
    #tprof_start = time_profile.perf_counter()

    # Temporary upscaled unmasked grid file.
    upscaled_unmasked_grd_file = tempfile.NamedTemporaryFile(delete=True)
    upscaled_unmasked_grd_file.close()  # cannot open twice on Windows - close before opening again

    try:
        #tprof_upscaled_unmasked_grid_start = time_profile.perf_counter()

        # Write the unmasked grid at the upscaled grid spacing.
        call_system_command([
                "gmt",
                "nearneighbor",
                xyz_filename,
                "-N8+m1", # Only require a single sample.
                "-S{}d".format(2 * grid_spacing),
                "-I{}".format(upscaled_grid_spacing),
                # Use GMT gridline registration since our input point grid has data points on the grid lines.
                # Gridline registration is the default so we don't need to force pixel registration...
                #"-r", # Force pixel registration since data points are at centre of cells.
                "-R{}/{}/{}/{}".format(-180, 180, -90, 90),
                "-G{}".format(upscaled_unmasked_grd_file.name)])
        
        #tprof_upscaled_unmasked_grid_end = time_profile.perf_counter()
        #print(f"  nearneighbor: {tprof_upscaled_unmasked_grid_end - tprof_upscaled_unmasked_grid_start:.2f} seconds")
        #tprof_upscaled_masked_grid_start = time_profile.perf_counter()

        # Since we used 'nearneighbour' above, it expanded outside the age grid mask.
        # So here we remove any grid values outside the age grid mask (using "OR" operator only looks at whether age grid has NaN values or not).
        call_system_command(["gmt", "grdmath", "-fg", upscaled_unmasked_grd_file.name, upscaled_mask_grd_filename, "OR", "=", grd_filename])

        #tprof_upscaled_masked_grid_end = time_profile.perf_counter()
        #print(f"  grdmath: {tprof_upscaled_masked_grid_end - tprof_upscaled_masked_grid_start:.2f} seconds")
    
    finally:
        os.unlink(upscaled_unmasked_grd_file.name)  # remove temp file (because we set 'delete=False')
    
    #tprof_end = time_profile.perf_counter()
    #print(f"write_upscaled_grd_file_from_xyz: {tprof_end - tprof_start:.2f} seconds")


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
            # Statistics to calculate mean/std-dev for a single ocean basin point.
            # Each ocean basin point has (num_proximities, sum_proximities, sum_square_proximities) that start at zero.
            self.point_statistics = np.full((len(self.point_lons), 3), (0, 0.0, 0.0), dtype=float)  # numpy array uses less memory
        
        # Keep a record of all proximity over time.
        if self.output_proximity_with_time:
            self.time_datas = {}  # dict indexed by time
    
    def add_proximity(self, proximity_in_kms, time, ocean_basin_point_index, ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat):
        # Update the proximity statistics for the current ocean basin point.
        if self.output_mean_proximity or self.output_standard_deviation_proximity:
            num_proximities, sum_proximities, sum_square_proximities = self.point_statistics[ocean_basin_point_index]
            num_proximities += 1
            sum_proximities += proximity_in_kms
            sum_square_proximities += proximity_in_kms * proximity_in_kms
            self.point_statistics[ocean_basin_point_index] = (num_proximities, sum_proximities, sum_square_proximities)
        
        # Add proximity for the current reconstructed point to a list for the reconstruction time.
        if self.output_proximity_with_time:
            # If we haven't already, create a new array of (reconstructed_lon, reconstructed_lat, proximity_in_kms) for all points at the specified time.
            if time not in self.time_datas:
                self.time_datas[time] = np.full((len(self.point_lons), 3), (0.0, 0.0, 0.0), dtype=float)  # numpy array uses less memory
            # Store the data for the current ocean point.
            self.time_datas[time][ocean_basin_point_index] = (ocean_basin_reconstructed_lon, ocean_basin_reconstructed_lat, proximity_in_kms)
    
    # Return the list of times (added with 'add_proximity()').
    def get_times(self):
        return list(self.time_datas.keys())
    
    # Return array of proximity tuples for the specified time.
    # Actually each tuple is an array of length 3 containing (reconstructed_lon, reconstructed_lat, proximity_in_kms).
    def get_time_data(self, time):
        if self.output_proximity_with_time:
            return self.time_datas[time]  # raises KeyError if time not in dict
        else:
            return []
    
    # Return list of mean statistics (over time) tuples.
    # Each tuple is (lon, lat, mean).
    def get_mean_data(self):
        if self.output_mean_proximity:
            # Calculate a mean proximity over time for each ocean basin point.
            mean_data = []
            for point_index, (num_proximities, sum_proximities, sum_square_proximities) in enumerate(self.point_statistics):
                # Only add current point if it didn't get deactivated immediately (and hence has no statistics).
                if num_proximities > 0:
                    mean_proximity = sum_proximities / num_proximities

                    # Clamp mean proximity if requested.
                    if (self.clamp_mean_proximity_in_kms is not None and
                        mean_proximity > self.clamp_mean_proximity_in_kms):
                        mean_proximity = self.clamp_mean_proximity_in_kms

                    ocean_basin_lon = self.point_lons[point_index]
                    ocean_basin_lat = self.point_lats[point_index]

                    mean_data.append((ocean_basin_lon, ocean_basin_lat, mean_proximity))

            return mean_data
        else:
            return []
    
    # Return list of standard deviation statistics (over time) tuples.
    # Each tuple is (lon, lat, std_dev).
    def get_std_dev_data(self):
        if self.output_standard_deviation_proximity:
            # Calculate a standard deviation proximity over time for each ocean basin point.
            std_dev_data = []
            for point_index, (num_proximities, sum_proximities, sum_square_proximities) in enumerate(self.point_statistics):
                # Only add current point if it didn't get deactivated immediately (and hence has no statistics).
                if num_proximities > 0:
                    mean_proximity = sum_proximities / num_proximities
                
                    standard_deviation_proximity_squared = (sum_square_proximities / num_proximities) - (mean_proximity * mean_proximity)
                    # Ensure not negative due to numerical precision.
                    if standard_deviation_proximity_squared > 0:
                        standard_deviation_proximity = math.sqrt(standard_deviation_proximity_squared)
                    else:
                        standard_deviation_proximity = 0
                        
                        ocean_basin_lon = self.point_lons[point_index]
                        ocean_basin_lat = self.point_lats[point_index]
                    
                    std_dev_data.append((ocean_basin_lon, ocean_basin_lat, standard_deviation_proximity))

            return std_dev_data
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
        shortest_path_grid = shortest_path.Grid(6)
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
                
                # If the age grid paleo time does not coincide with a time step then reconstruct from the former to the latter.
                # This can happen if the time interval is larger than 1 Myr.
                if time > age_grid_paleo_time:
                    # Resolve our topological plate polygons at the age grid paleo time.
                    resolved_plate_boundaries = []
                    pygplates.resolve_topologies(
                            topology_reconstruction_features,
                            rotation_model,
                            resolved_plate_boundaries,
                            age_grid_paleo_time,
                            resolve_topology_types=pygplates.ResolveTopologyType.boundary)
                    resolved_plate_polygons = [resolved_plate_boundary.get_resolved_boundary()
                            for resolved_plate_boundary in resolved_plate_boundaries]
                    # Reconstruct the current ocean basin points from 'age_grid_paleo_time' to 'time'.
                    topology_reconstruct_time_step(
                            ocean_basin_reconstruction,
                            age_grid_paleo_time,
                            time - age_grid_paleo_time,  # time increment
                            resolved_plate_boundaries,
                            resolved_plate_polygons,
                            rotation_model)
                    del resolved_plate_boundaries  # free memory
                    #print('Reconstructed initial age grid {} to time {}'.format(age_grid_paleo_time, time))
                
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
            topology_obstacle_resolved_topologies = []
            topology_obstacle_shared_boundary_sections = []
            pygplates.resolve_topologies(topology_reconstruction_features, rotation_model, topology_obstacle_resolved_topologies, time, topology_obstacle_shared_boundary_sections)
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
            del topology_obstacle_resolved_topologies
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

            cpu_profile.start_topology_resolve_time_step()
            
            # Resolve our topological plate polygons.
            resolved_plate_boundaries = []
            pygplates.resolve_topologies(
                    topology_reconstruction_features,
                    rotation_model,
                    resolved_plate_boundaries,
                    time,
                    resolve_topology_types=pygplates.ResolveTopologyType.boundary)
            resolved_plate_polygons = [resolved_plate_boundary.get_resolved_boundary()
                    for resolved_plate_boundary in resolved_plate_boundaries]
    
            cpu_profile.end_topology_resolve_time_step()
            
            # Find the ocean basin points for the next time step.
            # We do this for each age grid currently being reconstructed.
            for age_grid_paleo_time in list(ocean_basin_reconstructions.keys()):  # copy dict keys since might remove them while iterating
                ocean_basin_reconstruction = ocean_basin_reconstructions[age_grid_paleo_time]
                # Reconstruct the current ocean basin points from 'time' to 'time + time_increment'.
                # The reconstructed points will be the current points in the next time step.
                topology_reconstruct_time_step(
                        ocean_basin_reconstruction,
                        time,
                        time_increment,
                        resolved_plate_boundaries,
                        resolved_plate_polygons,
                        rotation_model)
                # If finished reconstructing ocean basin (for associated age grid) then remove from current reconstructions.
                if not ocean_basin_reconstruction.is_active():
                    #print('Finished age grid {} at time {}'.    format(age_grid_paleo_time, time))
                    del ocean_basin_reconstructions[age_grid_paleo_time]
    
            del resolved_plate_boundaries  # free memory
    
        cpu_profile.end_reconstruct_time_step()
        
        # Increment the time (to the next time interval).
        time_index += 1

    cpu_profile.end_reconstruct_and_calculate_distances()
    cpu_profile.end_proximity()
    cpu_profile.print_proximity_usage(age_grid_paleo_times)
    
    #memory_profile.print_object_memory_usage(proximity_datas, 'proximity_datas')
    for proximity_data_time in proximity_datas.keys():
        memory_profile.print_object_memory_usage(proximity_datas[proximity_data_time], 'proximity_datas[{}]'.format(proximity_data_time))

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
    
    
    # Write the distance grid(s) associated with each input age grid.
    for age_grid_filename, age_grid_paleo_time in age_grid_filenames_and_paleo_times:

        proximity_data = proximity_datas.get(age_grid_paleo_time)  # lookup in dict
        if not proximity_data:
            print('WARNING: All ocean basin points are outside the age grid: {}'.format(age_grid_filename), file=sys.stderr)
            continue

        # Create a temporary mask grid file (matching age grid mask) at our upscaled grid spacing.
        # We only do this if we're outputting mean/std-dev grids, and if upscaling has been enabled.
        upscaled_mask_grd_file = None
        if output_mean_distance or output_standard_deviation_distance:
            if output_grd_files:
                _, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    upscaled_mask_grd_file = get_upscaled_mask_grd_file(upscale_mean_std_dev_grid_spacing, age_grid_filename)
    
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
        
        if output_mean_distance:

            xyz_mean_distance_filename = '{}_{:.1f}_mean_distance.{}'.format(output_filename_prefix, age_grid_paleo_time, output_filename_extension)
            write_xyz_file(xyz_mean_distance_filename, proximity_data.get_mean_data())
            
            if output_grd_files:
                grd_mean_distance_filename = '{}_{:.1f}_mean_distance.nc'.format(output_filename_prefix, age_grid_paleo_time)
                ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    write_upscaled_grd_file_from_xyz(
                            grd_mean_distance_filename, xyz_mean_distance_filename,
                            ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing,
                            upscaled_mask_grd_file.name)
                else:
                    write_grd_file_from_xyz(
                            grd_mean_distance_filename, xyz_mean_distance_filename, ocean_basin_grid_spacing,
                            # Using original (grid-aligned) points so don't near nearest neighbour filtering...
                            use_nearneighbor=False)
        
        if output_standard_deviation_distance:

            xyz_standard_deviation_distance_filename = '{}_{:.1f}_std_dev_distance.{}'.format(output_filename_prefix, age_grid_paleo_time, output_filename_extension)
            write_xyz_file(xyz_standard_deviation_distance_filename, proximity_data.get_std_dev_data())

            if output_grd_files:
                grd_standard_deviation_distance_filename = '{}_{:.1f}_std_dev_distance.nc'.format(output_filename_prefix, age_grid_paleo_time)
                ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing = output_grd_files
                if upscale_mean_std_dev_grid_spacing is not None:
                    write_upscaled_grd_file_from_xyz(
                            grd_standard_deviation_distance_filename, xyz_standard_deviation_distance_filename,
                            ocean_basin_grid_spacing, upscale_mean_std_dev_grid_spacing,
                            upscaled_mask_grd_file.name)
                else:
                    write_grd_file_from_xyz(
                            grd_standard_deviation_distance_filename, xyz_standard_deviation_distance_filename, ocean_basin_grid_spacing,
                            # Using original (grid-aligned) points so don't near nearest neighbour filtering...
                            use_nearneighbor=False)
        
        # Remove temporary mask grid file (if we created it).
        if upscaled_mask_grd_file is not None:
            os.unlink(upscaled_mask_grd_file.name)  # remove temp file (because we set 'delete=False')
    
    # See how much extra memory is used after time/mean/std-dev data is extracted from the ProximityData.
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

    write_proximity_data(
            proximity_datas,
            age_grid_filenames_and_paleo_times,
            output_filename_prefix,
            output_filename_extension,
            output_distance_with_time,
            output_mean_distance,
            output_standard_deviation_distance,
            output_grd_files)


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
    base_memory_usage_per_task_in_gb = 2.0
    # The memory usage per age grid is roughly proportional to the number of input points,
    # with a uniform lon-lat grid at 1 degree resolution consuming about 6MB.
    delta_memory_usage_per_age_grid_in_gb = 6e-3 * len(input_points) / (180 * 360)
    # The total memory used to process the specified number of age grids in a single task.
    def memory_usage_per_task(num_age_grids_per_task_):
        return base_memory_usage_per_task_in_gb + num_age_grids_per_task_ * delta_memory_usage_per_age_grid_in_gb

    # If we've been given a limit on memory usage then determine how many age grids to process per task.
    if max_memory_usage_in_gb:
        # The memory used by the number of age grids per task multiplied by the number of tasks processed in parallel should not exceed the maximum memory usage.
        num_age_grids_per_task = math.trunc(((max_memory_usage_in_gb / num_cpus) - base_memory_usage_per_task_in_gb) / delta_memory_usage_per_age_grid_in_gb)  # could be negative
        # But don't reduce below the minimum number of age grids per task.
        if num_age_grids_per_task < min_num_age_grids_per_task:
            num_age_grids_per_task = min_num_age_grids_per_task
            # Reduce the number of CPUs to compensate for the higher than expected number of age grids per task, so that we don't exceed the memory limit.
            # Number of CPUs is the max memory divided by the memory used to process 'num_age_grids_per_task' age grids.
            num_cpus = math.trunc(max_memory_usage_in_gb / memory_usage_per_task(num_age_grids_per_task))
            if num_cpus < 1:
                num_cpus = 1
    else:
        # No limits on memory usage were specified, so just use the minimum number of age grids per task.
        num_age_grids_per_task = min_num_age_grids_per_task
    
    # If there are fewer total tasks than twice the number of CPUs then reduce the number of age grids per task
    # (but not less than the minimum) so that each CPU gets two tasks to process.
    # This helps to better utilise all CPUs when each task takes a different amount of time to complete.
    num_total_tasks = math.ceil(num_age_grids / num_age_grids_per_task)
    if num_total_tasks < 2 * num_cpus:
        # Note that reducing the number of age grids per task uses less memory (ie, doesn't violate our max memory limit).
        num_age_grids_per_task = math.ceil(num_age_grids / (2 * num_cpus))
        if num_age_grids_per_task < min_num_age_grids_per_task:
            num_age_grids_per_task = min_num_age_grids_per_task
        num_total_tasks = math.ceil(num_age_grids / num_age_grids_per_task)
        # If there are fewer tasks than the number of CPUs reduce the number of CPUs to match.
        if num_total_tasks < num_cpus:
            num_cpus = num_total_tasks

    if max_memory_usage_in_gb:
        print('Maximum memory usage: {:.2f}GB'.format(max_memory_usage_in_gb))
    print('Approximate memory usage: {:.2f}GB'.format(num_cpus * memory_usage_per_task(num_age_grids_per_task)))
    print('Number of age grids: {}'.format(num_age_grids))
    print('Number of age grids per task: {}'.format(num_age_grids_per_task))
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
