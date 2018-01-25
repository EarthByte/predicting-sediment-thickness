
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


from __future__ import print_function
import argparse
from call_system_command import call_system_command
import math
import multiprocessing
import points_in_polygons
import proximity_query
import pygplates
import subprocess
import sys


# Default grid spacing (in degrees) when generating uniform lon/lat spacing of ocean basin points.
DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES = 1


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
# For input points outside the age grid then ages will be Nan (ie, 'math.isnan(age)' will return True).
def get_positions_and_ages(input_points, age_grid_filename):
    
    input_points_data = ''.join('{0} {1}\n'.format(lon, lat) for lon, lat in input_points)
    
    stdout_data = call_system_command(
            # The command-line strings to execute GMT 'grdtrack'...
            ["gmt", "grdtrack", "-nl", "-G{0}".format(age_grid_filename)],
            stdin=input_points_data,
            return_stdout=True)
    
    #print('Stdout: {0}'.format(stdout_data))
    
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
            print('Ignoring line "{0}" - has fewer than 3 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            lon = float(line_data[0])
            lat = float(line_data[1])
            
            # The age got appended to the last column by 'grdtrack'.
            age = float(line_data[-1])
            
            # If the point is outside the ocean basin region then the age grid will return 'NaN'.
            if math.isnan(age):
                #print('Ignoring line "{0}" - point is outside ocean basin (age grid).'.format(line), file=sys.stderr)
                continue
            
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point lon, lat and age values.'.format(line), file=sys.stderr)
            continue
        
        lon_lat_age_list.append((lon, lat, age))
    
    return lon_lat_age_list


# Determine the overriding plate of the subducting line.
def find_overriding_plate(subduction_shared_sub_segment, time, reconstructed_point):
    
    # Get the subduction polarity of the nearest subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
        print('Unable to find the overriding plate of the nearest subducting line "{0}"'.format(
            subduction_shared_sub_segment.get_feature().get_name()), file=sys.stderr)
        print('    subduction zone feature is missing subduction polarity property or it is set to "Unknown".', file=sys.stderr)
        return

    overriding_plate = None

    # Iterate over the topologies that are sharing the part (sub-segment) of the subducting line that is closest to the feature.
    sharing_resolved_topologies = subduction_shared_sub_segment.get_sharing_resolved_topologies()
    geometry_reversal_flags = subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    for index in range(len(sharing_resolved_topologies)):

        sharing_resolved_topology = sharing_resolved_topologies[index]
        geometry_reversal_flag = geometry_reversal_flags[index]

        if sharing_resolved_topology.get_resolved_boundary().get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise:
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if ((subduction_polarity == 'Left' and geometry_reversal_flag) or
                (subduction_polarity == 'Right' and not geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
                break
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is not reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if ((subduction_polarity == 'Left' and not geometry_reversal_flag) or
                (subduction_polarity == 'Right' and geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
                break
    
    if overriding_plate is None:
        print('Unable to find the overriding plate of the nearest subducting line "{0}" at {1}Ma for reconstructed point {2}'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time, reconstructed_point.to_lat_lon()), file=sys.stderr)
        print('    topology on overriding side of subducting line is missing.', file=sys.stderr)
        return
    
    return overriding_plate


# Reconstruct points from 'time' to 'time + time_increment' using topologies.
def topology_reconstruct_time_step(
        time,
        time_increment,
        resolved_plate_boundaries,
        stage_rotation_cache,
        curr_point_infos,
        next_point_infos):
    #
    # The following code is the *slower* version of point-in-polygon testing.
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
    #        stage_rotation = stage_rotation_cache.get_stage_rotation(plate_id, time, time_increment)
    #        next_recon_point = stage_rotation * curr_recon_point
    #        
    #        next_point_infos.append((point_feature, point_begin_time, next_recon_point))
    
    #
    # The following code is the *faster* version of the above point-in-polygon code.
    #
    
    # Remove any points that don't exist at 'time + time_increment'.
    curr_existing_point_infos = []
    curr_existing_recon_points = []
    for curr_point_info in curr_point_infos:
        point_begin_time = curr_point_info[1]
        # Retire current point if the time we are reconstructing to ('time + time_increment')
        # is older (earlier than) than the point's time of appearance.
        if time + time_increment > point_begin_time:
            continue
        
        curr_existing_recon_point = curr_point_info[2]
        curr_existing_recon_points.append(curr_existing_recon_point)
        
        curr_existing_point_infos.append(curr_point_info)
        
    # Find the polygon plates containing the points.
    resolved_plate_polygons = [resolved_plate_boundary.get_resolved_boundary()
            for resolved_plate_boundary in resolved_plate_boundaries]
    resolved_plate_boundaries_containing_points = points_in_polygons.find_polygons(
            curr_existing_recon_points, resolved_plate_polygons, resolved_plate_boundaries)
    
    # Iterate over the point-in-polygon results.
    for point_index, resolved_plate_boundary in enumerate(resolved_plate_boundaries_containing_points):
        if resolved_plate_boundary is None:
            continue
        
        curr_paleo_lon_lat_point, point_begin_time, curr_recon_point = curr_existing_point_infos[point_index]
        
        plate_id = resolved_plate_boundary.get_feature().get_reconstruction_plate_id()
        
        # Get the stage rotation from 'time' to 'time + time_increment'.
        stage_rotation = stage_rotation_cache.get_stage_rotation(plate_id, time, time_increment)
        next_recon_point = stage_rotation * curr_recon_point
        
        next_point_infos.append((curr_paleo_lon_lat_point, point_begin_time, next_recon_point))


def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def write_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, num_grid_longitudes, num_grid_latitudes):
    
    # The command-line strings to execute GMT 'nearneighbor'.
    # For example "nearneighbor output_mean_distance.xy -R-179.5/179.5/-89.5/89.5 -I1 -N4 -S1d -Goutput_mean_distance.grd=cf"
    # ...with "=cf" generating NetCDF-3, instead of NetCDF-4, files (GPlates 2.0 can only load NetCDF-3).
    gmt_command_line = [
            "gmt",
            "nearneighbor",
            xyz_filename.encode(sys.getfilesystemencoding()),
            "-N4",
            "-S{0}d".format(1.5 * grid_spacing),
            "-I{0}".format(grid_spacing),
            # Use GMT gridline registration since our input point grid has data points on the grid lines.
            # Gridline registration is the default so we don't need to force pixel registration...
            #"-r", # Force pixel registration since data points are at centre of cells.
            "-R{0}/{1}/{2}/{3}".format(-180, 180, -90, 90),
            #"-R{0}/{1}/{2}/{3}".format(
            #        -180 + 0.5 * grid_spacing,
            #        -180 + (num_grid_longitudes - 0.5) * grid_spacing,
            #        -90 + 0.5 * grid_spacing,
            #        -90 + (num_grid_latitudes - 0.5) * grid_spacing),
            "-G{0}=cf".format(grd_filename.encode(sys.getfilesystemencoding()))]
    
    call_system_command(gmt_command_line)


# Class to calculate and re-use stage rotation calculations.
class StageRotationCache(object):
    def __init__(self, rotation_model, anchor_plate_id):
        self.rotation_model = rotation_model
        self.anchor_plate_id = anchor_plate_id
        self.stage_rotation_dict = {}
    
    def get_stage_rotation(self, plate_id, time, time_increment):
        stage_rotation = self.stage_rotation_dict.get((plate_id, time, time_increment))
        if not stage_rotation:
            stage_rotation = self.rotation_model.get_rotation(
                    time + time_increment, plate_id, time, anchor_plate_id=self.anchor_plate_id)
            # Cache for next time.
            self.stage_rotation_dict[(plate_id, time, time_increment)] = stage_rotation
        
        return stage_rotation


# Class to hold all proximity data associated with a particular proximity feature.
class FeatureProximityData(object):
    def __init__(self):
        self.time_data = {}
        self.mean_data = []
        self.std_dev_data = []
    
    # Return list of data tuples at specified time.
    # Each tuple is (lon, lat [, proximity] [, overriding_plate_id])
    # to write to the output file for the current 'time'.
    def get_time_data(self, time_index):
        return self.time_data.setdefault(time_index, [])
    
    # Return dict of all time data.
    def get_all_time_data(self):
        return self.time_data
    
    # Return list of mean statistics (over time) tuples.
    # Each tuple is (lon, lat, mean)
    def get_mean_data(self):
        return self.mean_data
    
    # Return list of standard deviation statistics (over time) tuples.
    # Each tuple is (lon, lat, std_dev)
    def get_std_dev_data(self):
        return self.std_dev_data
    
    # Merge other feature proximity data in us.
    def update(self, other):
        for time_index, other_time_data in other.get_all_time_data().iteritems():
            time_data = self.get_time_data(time_index)
            time_data += other_time_data
        self.mean_data += other.get_mean_data();
        self.std_dev_data += other.get_std_dev_data();


# Class to hold all proximity data.
class ProximityData(object):
    def __init__(self):
        self.feature_data = {}
    
    # Return FeatureProximityData object for specified feature name.
    def get_feature_data(self, feature_name=None):
        return self.feature_data.setdefault(feature_name, FeatureProximityData())
    
    # Return dict of all feature data.
    def get_all_feature_data(self):
        return self.feature_data
    
    # Merge other proximity data into us.
    def update(self, other):
        for feature_name, other_feature_data in other.get_all_feature_data().iteritems():
            self.get_feature_data(feature_name).update(other_feature_data)


def proximity_parallel(
    input_points, # List of (lon, lat) tuples.
    rotation_filenames,
    proximity_filenames,
    proximity_features_are_topological,
    proximity_feature_types,
    topological_reconstruction_filenames,
    max_topological_reconstruction_time,
    age_grid_filename,
    age_grid_paleo_time,
    time_increment,
    output_distance_with_time,
    output_overriding_plate_with_time,
    output_mean_distance,
    output_standard_deviation_distance,
    output_all_proximity_distances,
    anchor_plate_id = 0,
    proximity_distance_threshold_radians = None):
    
    if time_increment <= 0:
        print('The time increment "{0}" is not positive and non-zero.'.format(time_increment), file=sys.stderr)
        return
    
    if age_grid_paleo_time < 0:
        print('The age grid paleo time "{0}" is not positive or zero.'.format(age_grid_paleo_time), file=sys.stderr)
        return
    
    if (not output_distance_with_time and
        not output_overriding_plate_with_time and
        not output_mean_distance and
        not output_standard_deviation_distance):
        print('No output specified.', file=sys.stderr)
        return
    
    if (output_all_proximity_distances and
        proximity_features_are_topological):
        print('Outputting all proximity distances not supported for topological features.', file=sys.stderr)
        return
    
    # Make sure user provides *topological* proximity features and specifies subduction zone when
    # they want to output overriding plate IDs.
    if output_overriding_plate_with_time:
        if not proximity_features_are_topological:
            return
        if (not proximity_feature_types or
            len(proximity_feature_types) != 1 or
            proximity_feature_types[0] != pygplates.FeatureType.gpml_subduction_zone.get_name()):
            return
    
    # Get the input point ages.
    lon_lat_age_list = get_positions_and_ages(input_points, age_grid_filename)
    if not lon_lat_age_list:
        return
    
    # For each ocean basin point (and associated age) create a tuple containing
    # ((paleo_lon, paleo_lat), time_of_appearance, paleo_point).
    ocean_basin_point_infos = []
    for lon, lat, age in lon_lat_age_list:
        # List of 3-tuples (point_feature, point_begin_time, point).
        ocean_basin_point_infos.append(
                ((lon, lat), age_grid_paleo_time + age, pygplates.PointOnSphere(lat, lon)))
    
    rotation_model = pygplates.RotationModel(rotation_filenames)
    
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
    
    stage_rotation_cache = StageRotationCache(rotation_model, anchor_plate_id)
    
    if output_mean_distance or output_standard_deviation_distance:
        # Create a dictionary mapping each ocean basin point (lat/lon position) to its proximity statistics.
        ocean_basin_points_statistics = { }
    
    # All proximity data to return to caller.
    proximity_data = ProximityData()
    
    # Iterate from paleo time until we exceed the maximum begin time of all ocean basin point locations.
    min_time_index = int(math.ceil(age_grid_paleo_time / time_increment))
    max_time = max(ocean_basin_point_info[1] for ocean_basin_point_info in ocean_basin_point_infos)
    # We cannot reconstruct further in the past than allowed by the topological reconstruction features.
    if max_time > max_topological_reconstruction_time:
        max_time = max_topological_reconstruction_time
    max_time_index = int(math.trunc(max_time / time_increment))
    next_ocean_basin_point_infos = ocean_basin_point_infos
    for time_index in range(min_time_index, max_time_index + 1):
        
        time = time_index * time_increment
        #print('Time {0}'.format(time))
            
        current_ocean_basin_point_infos = next_ocean_basin_point_infos
        
        # Reconstruct to the next time unless we're already at the last time.
        if time_index != max_time_index:
            # Resolve our topological plate polygons.
            resolved_plate_boundaries = []
            pygplates.resolve_topologies(
                    topology_reconstruction_features,
                    rotation_model,
                    resolved_plate_boundaries,
                    time,
                    resolve_topology_types=pygplates.ResolveTopologyType.boundary)
            
            next_ocean_basin_point_infos = []
            
            # Reconstruct the current ocean basin points from 'time' to 'time + time_increment'.
            # The reconstructed points will be the current points in the next time step.
            topology_reconstruct_time_step(
                    time,
                    time_increment,
                    resolved_plate_boundaries,
                    stage_rotation_cache,
                    current_ocean_basin_point_infos,
                    next_ocean_basin_point_infos)
        
        ocean_basin_reconstructed_points = [current_ocean_basin_point_info[2]
                for current_ocean_basin_point_info in current_ocean_basin_point_infos]
        
        if proximity_features_are_topological:
            # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
            # We generate both the resolved topology boundaries and the boundary sections between them.
            proximity_resolved_topologies = []
            proximity_shared_boundary_sections = []
            pygplates.resolve_topologies(proximity_features, rotation_model, proximity_resolved_topologies, time, proximity_shared_boundary_sections, anchor_plate_id)
            
            # 'output_all_proximity_distances' is required to be False for topological proximity features.
            proximity_reconstructed_geometry_proxies = None
            
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
            
        else: # non-topological features...
            
            # Reconstruct the non-topological features that exist at the current 'time'.
            proximity_reconstructed_feature_geometries = []
            pygplates.reconstruct(proximity_features, rotation_model, proximity_reconstructed_feature_geometries, time, anchor_plate_id)
            
            if output_all_proximity_distances:
                proximity_reconstructed_geometry_proxies = []
            else:
                proximity_reconstructed_geometry_proxies = None
            
            proximity_reconstructed_geometries = []
            for proximity_reconstructed_feature_geometry in proximity_reconstructed_feature_geometries:
                proximity_reconstructed_geometries.append(proximity_reconstructed_feature_geometry.get_reconstructed_geometry())
                if proximity_reconstructed_geometry_proxies is not None:
                    proximity_reconstructed_geometry_proxies.append(proximity_reconstructed_feature_geometry.get_feature())
        
        # Find the minimum distance of each ocean basin point to all proximity reconstructed geometries.
        proximity_features_closest_to_ocean_basin_points = proximity_query.find_closest_geometries_to_points(
                ocean_basin_reconstructed_points,
                proximity_reconstructed_geometries,
                proximity_reconstructed_geometry_proxies,
                distance_threshold_radians = proximity_distance_threshold_radians,
                # Does each point get a list of all geometries within threshold distance (instead of just the closest) ? ...
                all_geometries = output_all_proximity_distances)
        
        # Iterate over all reconstructed ocean basin points.
        for ocean_basin_point_index, proximity_feature_closest_to_ocean_basin_point in enumerate(proximity_features_closest_to_ocean_basin_points):
            if proximity_feature_closest_to_ocean_basin_point is None:
                # All geometries were further than the distance threshold, so skip current point.
                continue
            
            # Mapping of proximity feature names to proximity distances (in radians).
            # The special key 'None' is associated with the *nearest* proximity distance over *all* proximity features.
            ocean_basin_proximities = {}
            
            if output_all_proximity_distances:
                # 'proximity_feature_closest_to_ocean_basin_point' is actually a list of all proximity features
                # within the distance threshold of the current ocean basin point (not just the closest proximity feature).
                proximity_features_closest_to_ocean_basin_point = proximity_feature_closest_to_ocean_basin_point
                
                # The minimum distance of the current ocean basin point to all reconstructed proximity features.
                min_distance_to_all_proximity_features = None
                for min_distance_to_proximity_feature, proximity_feature in proximity_features_closest_to_ocean_basin_point:
                    proximity_feature_name = proximity_feature.get_name()
                    # Set proximity to the minimum distance to all features with the same name (if any).
                    # A proximity feature may have more than one geometry (we choose the closest one for that feature).
                    if proximity_feature_name in ocean_basin_proximities:
                        if min_distance_to_proximity_feature < ocean_basin_proximities[proximity_feature_name]:
                            ocean_basin_proximities[proximity_feature_name] = min_distance_to_proximity_feature
                    else:
                        ocean_basin_proximities[proximity_feature_name] = min_distance_to_proximity_feature
                    
                    # If the current feature is the closest of all the features so far.
                    if (min_distance_to_all_proximity_features is None or
                            min_distance_to_proximity_feature < min_distance_to_all_proximity_features):
                        min_distance_to_all_proximity_features = min_distance_to_proximity_feature
                
                # Minimum distance to *all* proximity features is represented by 'None'.
                ocean_basin_proximities[None] = min_distance_to_all_proximity_features
            else:
                min_distance_to_all_proximity_features, _ = proximity_feature_closest_to_ocean_basin_point
                # Minimum distance to *all* proximity features is represented by 'None'.
                ocean_basin_proximities[None] = min_distance_to_all_proximity_features
            
            # Extract the current ocean basin point parameters.
            ocean_basin_paleo_lon_lat_point, ocean_basin_point_begin_time, ocean_basin_reconstructed_point = current_ocean_basin_point_infos[ocean_basin_point_index]
            
            if output_mean_distance or output_standard_deviation_distance:
                # Update the distance statistics for the current ocean basin point.
                ocean_basin_point_statistics = ocean_basin_points_statistics.setdefault(ocean_basin_paleo_lon_lat_point, {})
                for proximity_feature_name, proximity in ocean_basin_proximities.iteritems():
                    proximity_in_kms = proximity * pygplates.Earth.mean_radius_in_kms
                    num_distances, sum_distances, sum_square_distances = ocean_basin_point_statistics.setdefault(proximity_feature_name, (0, 0.0, 0.0))
                    num_distances += 1
                    sum_distances += proximity_in_kms
                    sum_square_distances += proximity_in_kms * proximity_in_kms
                    ocean_basin_point_statistics[proximity_feature_name] = (num_distances, sum_distances, sum_square_distances)
            
            if output_distance_with_time or output_overriding_plate_with_time:
                # Determine the overriding plate of the subducting line (if requested).
                if output_overriding_plate_with_time:
                    overriding_plate = find_overriding_plate(
                            nearest_proximity_shared_sub_segment,
                            time,
                            ocean_basin_reconstructed_point)
                    if not overriding_plate:
                        # The overriding plate was not found so skip the current ocean basin point.
                        # Note that a warning has already been generated by "find_overriding_plate()".
                        continue
                
                ocean_basin_reconstructed_lat, ocean_basin_reconstructed_lon = ocean_basin_reconstructed_point.to_lat_lon()
                
                for proximity_feature_name, proximity in ocean_basin_proximities.iteritems():
                    # Group first by feature name and then by time.
                    # This makes it easier to do the final writes to the output files.
                    feature_proximity_data = proximity_data.get_feature_data(proximity_feature_name)
                    feature_proximity_time_data = feature_proximity_data.get_time_data(time_index)
                    
                    if output_distance_with_time:
                        if output_overriding_plate_with_time:
                            feature_proximity_time_data.append((
                                ocean_basin_reconstructed_lon,
                                ocean_basin_reconstructed_lat,
                                proximity * pygplates.Earth.mean_radius_in_kms,
                                overriding_plate.get_feature().get_reconstruction_plate_id()))
                        else:
                            feature_proximity_time_data.append((
                                ocean_basin_reconstructed_lon,
                                ocean_basin_reconstructed_lat,
                                proximity * pygplates.Earth.mean_radius_in_kms))
                    else:
                        feature_proximity_time_data.append((
                            ocean_basin_reconstructed_lon,
                            ocean_basin_reconstructed_lat,
                            overriding_plate.get_feature().get_reconstruction_plate_id()))
    
    if output_mean_distance or output_standard_deviation_distance:
        for (ocean_basin_lon, ocean_basin_lat), ocean_basin_point_statistics in ocean_basin_points_statistics.iteritems():
            
            for proximity_feature_name, (num_distances, sum_distances, sum_square_distances) in ocean_basin_point_statistics.iteritems():
                # Group by feature name. This makes it easier to do the final writes to the output files.
                feature_proximity_data = proximity_data.get_feature_data(proximity_feature_name)
                
                mean_distance = sum_distances / num_distances
                
                if output_mean_distance:
                    feature_proximity_data.get_mean_data().append((ocean_basin_lon, ocean_basin_lat, mean_distance))
                
                if output_standard_deviation_distance:
                    standard_deviation_distance_squared = (sum_square_distances / num_distances) - (mean_distance * mean_distance)
                    # Ensure not negative due to numerical precision.
                    if standard_deviation_distance_squared > 0:
                        standard_deviation_distance = math.sqrt(standard_deviation_distance_squared)
                    else:
                        standard_deviation_distance = 0
                    
                    feature_proximity_data.get_std_dev_data().append((ocean_basin_lon, ocean_basin_lat, standard_deviation_distance))
    
    return proximity_data


# Wraps around 'proximity()' so can be used by multiprocessing.Pool.map() which
# requires a single-argument function.
def proximity_parallel_pool_function(args):
    try:
        return proximity_parallel(*args)
    except KeyboardInterrupt:
        pass


def proximity(
        input_points, # List of (lon, lat) tuples.
        rotation_filenames,
        proximity_filenames,
        proximity_features_are_topological,
        proximity_feature_types,
        topological_reconstruction_filenames,
        max_topological_reconstruction_time,
        age_grid_filename,
        age_grid_paleo_time,
        time_increment,
        output_distance_with_time,
        output_overriding_plate_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        output_all_proximity_distances,
        num_cpus, # If None then defaults to all available CPUs.
        anchor_plate_id = 0,
        proximity_distance_threshold_kms = None):
    
    if time_increment <= 0:
        print('The time increment "{0}" is not positive and non-zero.'.format(time_increment), file=sys.stderr)
        return
    
    if age_grid_paleo_time < 0:
        print('The age grid paleo time "{0}" is not positive or zero.'.format(age_grid_paleo_time), file=sys.stderr)
        return
    
    if (not output_distance_with_time and
        not output_overriding_plate_with_time and
        not output_mean_distance and
        not output_standard_deviation_distance):
        print('No output specified.', file=sys.stderr)
        return
    
    if (output_all_proximity_distances and
        proximity_features_are_topological):
        print('Outputting all proximity distances not supported for topological features.', file=sys.stderr)
        return
    
    # Make sure user provides *topological* proximity features and specifies subduction zone when
    # they want to output overriding plate IDs.
    if output_overriding_plate_with_time:
        if not proximity_features_are_topological:
            return
        if (len(proximity_feature_types) != 1 or
            proximity_feature_types[0] != pygplates.FeatureType.gpml_subduction_zone.get_name()):
            return
    
    # If the user requested all available CPUs then attempt to find out how many there are.
    if not num_cpus:
        try:
            num_cpus = multiprocessing.cpu_count()
        except NotImplementedError:
            num_cpus = 1
    
    # Convert maximum proximity distance from Kms to radians (if it was specified).
    if proximity_distance_threshold_kms is None:
        proximity_distance_threshold_radians = None
    else:
        proximity_distance_threshold_radians = proximity_distance_threshold_kms / pygplates.Earth.mean_radius_in_kms
        if proximity_distance_threshold_radians > 2 * math.pi:
            # Exceeds circumference of Earth so no need for threshold.
            proximity_distance_threshold_radians = None
    
    #
    # Can uncomment this when there are exceptions in order to determine which source code line.
    #
    # Once goes through multiprocessing pools then lose error locations in source code.
    #
    
    # return proximity_parallel(
        # input_points, # List of (lon, lat) tuples.
        # rotation_filenames,
        # proximity_filenames,
        # proximity_features_are_topological,
        # proximity_feature_types,
        # topological_reconstruction_filenames,
        # max_topological_reconstruction_time,
        # age_grid_filename,
        # age_grid_paleo_time,
        # time_increment,
        # output_distance_with_time,
        # output_overriding_plate_with_time,
        # output_mean_distance,
        # output_standard_deviation_distance,
        # output_all_proximity_distances,
        # anchor_plate_id,
        # proximity_distance_threshold_radians)
    
    # Divide the input points into sub-lists (each to be run on a separate process).
    # Give each task a reasonable number of points to process - if there's not enough input points
    # per task then we'll spend too much time resolving/reconstructing proximity features
    # (which needs to be repeated for each task of input points).
    # So we could reduce the number of tasks (increase input points per task) to the number of CPUs.
    # However some tasks might finish sooner than others leaving some CPUs under utilised.
    # So we double the number of tasks (twice number of CPUs) to counteract this to some extent.
    num_tasks = max(2 * num_cpus, int(len(input_points) / 25000))
    #print('num_tasks', num_tasks)
    
    # Create the pool sub-lists.
    # The points in each sub-list are spread across the globe to help avoid some sub-lists
    # falling mostly outside ocean basin regions, and hence not getting processed, resulting
    # in some tasks progressing much faster than others.
    pool_input_points_sub_lists = [input_points[i:len(input_points):num_tasks]
            for i in xrange(num_tasks)]
    
    # Split the workload across the CPUs.
    pool = multiprocessing.Pool(num_cpus)
    pool_map_async_result = pool.map_async(
            proximity_parallel_pool_function,
            (
                (
                    pool_input_points_sub_list,
                    rotation_filenames,
                    proximity_filenames,
                    proximity_features_are_topological,
                    proximity_feature_types,
                    topological_reconstruction_filenames,
                    max_topological_reconstruction_time,
                    age_grid_filename,
                    age_grid_paleo_time,
                    time_increment,
                    output_distance_with_time,
                    output_overriding_plate_with_time,
                    output_mean_distance,
                    output_standard_deviation_distance,
                    output_all_proximity_distances,
                    anchor_plate_id,
                    proximity_distance_threshold_radians
                ) for pool_input_points_sub_list in pool_input_points_sub_lists
            ),
            1) # chunksize
    
    # Apparently if we use pool.map_async instead of pool.map and then get the results
    # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
    # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    try:
        pool_proximity_datas = pool_map_async_result.get(99999)
    except KeyboardInterrupt:
        return
    
    #
    # Extract and merge all proximity data from the pools.
    #
    
    proximity_data = ProximityData()
    
    for pool_proximity_data in pool_proximity_datas:
        if pool_proximity_data is None:
            return
        proximity_data.update(pool_proximity_data)
    
    return proximity_data;
    
    
def write_proximity_data(
        proximity_data,
        output_filename_prefix,
        output_filename_extension,
        time_increment,
        output_distance_with_time,
        output_overriding_plate_with_time,
        output_mean_distance,
        output_standard_deviation_distance,
        output_grd_files = None):
    
    all_proximity_feature_data = proximity_data.get_all_feature_data()
    for feature_name, proximity_feature_data in all_proximity_feature_data.iteritems():
        
        if output_distance_with_time or output_overriding_plate_with_time:
            for time_index, proximity_feature_time_data in proximity_feature_data.get_all_time_data().iteritems():
                time = time_index * time_increment
                if feature_name is not None:
                    xyz_filename = u'{0}_{1}_{2:0.2f}.{3}'.format(output_filename_prefix, feature_name.decode('utf-8'), time, output_filename_extension)
                    if output_grd_files:
                        grd_filename = u'{0}_{1}_{2:0.2f}.grd'.format(output_filename_prefix, feature_name.decode('utf-8'), time)
                else:
                    xyz_filename = u'{0}_{1:0.2f}.{2}'.format(output_filename_prefix, time, output_filename_extension)
                    if output_grd_files:
                        grd_filename = u'{0}_{1:0.2f}.grd'.format(output_filename_prefix, time)
                
                write_xyz_file(xyz_filename, proximity_feature_time_data)
                if output_grd_files:
                    grid_spacing, num_grid_longitudes, num_grid_latitudes = output_grd_files
                    write_grd_file_from_xyz(grd_filename, xyz_filename, grid_spacing, num_grid_longitudes, num_grid_latitudes)
        
        if output_mean_distance:
            if feature_name is not None:
                xyz_mean_distance_filename = u'{0}_mean_distance_{1}.{2}'.format(output_filename_prefix, feature_name.decode('utf-8'), output_filename_extension)
                if output_grd_files:
                    grd_mean_distance_filename = u'{0}_mean_distance_{1}.grd'.format(output_filename_prefix, feature_name.decode('utf-8'))
            else:
                xyz_mean_distance_filename = u'{0}_mean_distance.{1}'.format(output_filename_prefix, output_filename_extension)
                if output_grd_files:
                    grd_mean_distance_filename = u'{0}_mean_distance.grd'.format(output_filename_prefix)
                
            write_xyz_file(xyz_mean_distance_filename, proximity_feature_data.get_mean_data())
            if output_grd_files:
                grid_spacing, num_grid_longitudes, num_grid_latitudes = output_grd_files
                write_grd_file_from_xyz(grd_mean_distance_filename, xyz_mean_distance_filename, grid_spacing, num_grid_longitudes, num_grid_latitudes)
        
        if output_standard_deviation_distance:
            if feature_name is not None:
                xyz_standard_deviation_distance_filename = u'{0}_std_dev_distance_{1}.{2}'.format(output_filename_prefix, feature_name.decode('utf-8'), output_filename_extension)
                if output_grd_files:
                    grd_standard_deviation_distance_filename = u'{0}_std_dev_distance_{1}.grd'.format(output_filename_prefix, feature_name.decode('utf-8'))
            else:
                xyz_standard_deviation_distance_filename = u'{0}_std_dev_distance.{1}'.format(output_filename_prefix, output_filename_extension)
                if output_grd_files:
                    grd_standard_deviation_distance_filename = u'{0}_std_dev_distance.grd'.format(output_filename_prefix)
                
            write_xyz_file(xyz_standard_deviation_distance_filename, proximity_feature_data.get_std_dev_data())
            if output_grd_files:
                grid_spacing, num_grid_longitudes, num_grid_latitudes = output_grd_files
                write_grd_file_from_xyz(grd_standard_deviation_distance_filename, xyz_standard_deviation_distance_filename, grid_spacing, num_grid_longitudes, num_grid_latitudes)
    
    return 0 # Success


if __name__ == '__main__':
    
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
    
    The output can be any or all of the following:
     1) For each time from present day to the maximum age of all ocean basin points (determined by the age grid) an output xyzw file is
        generated containing the reconstructed (lon, lat) point locations (as x,y) and the minimum distance to all proximity features (as z) and
        optionally the overriding plate ID of the nearest subducting line (as w, or as z if '-d' option not specified) if the
        proximity features are topological and the proximity feature type ('-b') is subducting.
        See the '-d' and '-p' options.
     2) An output xyz file is generated containing 'present day' (lon, lat) ocean basin point locations (as x,y) and the mean distances
        of each point to all proximity features averaged over the point's lifetime (as determined by the age grid). See the '-j' option.
        A similar output xyz file can be generated containing standard deviation (instead of mean) distance. See the '-k' option.
    
    In addition to (1) and (2), equivalent xyz output files can be generated for each proximity feature. These are proximities to individual
    features (rather than proximity to the nearest of all proximity feature as in (1) and (2)). See the '-u' option.
    And a GMT grd file can be generated for each xyz file. See the '-w' option.
    
    If an ocean basin point falls outside the age grid then it is ignored. If the overriding plate of the nearest subduction line
    (if requested) cannot be found for a particular geological time then the point is ignored for that time step
    (in the output case (1) above). It is not ignored for output case (2) above (which does not involve overriding plates).
    
    A threshold distance can be specified to reject proximities exceeding it. See the '-q' option.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -m topologies.gpml -b SubductionZone -s static_polygons.gpml -g age_grid.nc -d -p -j -k -- input.xy output
     """
    
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
    parser.add_argument('-s', '--topological_reconstruction_filenames', type=str, nargs='+', required=True,
            metavar='topological_reconstruction_filename',
            help='The filenames of the topological files used to incrementally reconstruct the (paleo) ocean basin points.')
    parser.add_argument('-x', '--max_topological_reconstruction_time', type=int, required=True,
            metavar='max_topological_reconstruction_time',
            help='The maximum (largest) time that the topological reconstruction files can reconstruct back to (in Ma). Value must be an integer.')
    parser.add_argument('-g', '--age_grid_filename', type=str, required=True,
            metavar='age_grid_filename',
            help='The age grid filename used to find the begin age of each ocean basin point.')
    parser.add_argument('-y', '--age_grid_paleo_time', type=int, required=True,
            metavar='age_grid_paleo_time',
            help='The time of the age grid in Ma. Value must be an integer.')
    parser.add_argument('-t', '--time_increment', type=int, default=1,
            help='The time increment in My. Value must be an integer. Defaults to 1 My.')
    parser.add_argument('-q', '--max_distance', type=float,
            help='The maximum distance in Kms. '
                'If specified then distances between ocean basin points and proximity features exceeding this threshold will be ignored '
                '(ocean basin point will not get output or will not contribute to mean / standar deviation), otherwise all distances are included.')
    parser.add_argument('-c', '--num_cpus', type=int,
            help='The number of CPUs to use for calculations. Defaults to all available CPUs.')
    
    parser.add_argument('-d', '--output_distance_with_time', action='store_true',
            help='For each input point at each time during its lifetime write its distance to the nearest feature. '
                'If no output options (excluding "output_all_proximity_distances") are specified then this one is used.')
    parser.add_argument('-p', '--output_overriding_plate', action='store_true',
            dest='output_overriding_plate_with_time',
            help='For each input point at each time during its lifetime write the overriding plate ID of the nearest subduction zone. '
                'Write the overriding plate ID to the fourth (w) column of the output if also writing distance (-d/--output_distance_with_time), '
                'otherwise write to the third (z) column. '
                'This option can only be specified if using topological features (ie, -n is not specified) '
                'and subduction zones are the only feature type (ie, -b SubductionZone). '
                'By default it is not written.')
    parser.add_argument('-j', '--output_mean_distance', action='store_true',
            help='For each input point write its mean distance to features averaged over its lifetime. '
                'By default it is not written.')
    parser.add_argument('-k', '--output_std_dev_distance', action='store_true',
            help='For each input point write its standard deviation of distances to features averaged over its lifetime. '
                'By default it is not written.')
    parser.add_argument('-u', '--output_all_proximity_distances', action='store_true',
            help='Also output distances to each feature (an extra distance for each feature name). '
                'This outputs extra files with a feature name embedded in each filename. '
                'This is in addition to outputting distances to nearest features. '
                'By default only distances to the nearest features are written. '
                'Can only be specified if "non_topological_proximity_features" is also specified '
                '(ie, only applies to non-topological features).')
    parser.add_argument('-w', '--output_grd_files', action='store_true',
            help='Also generate a grd file for each xyz file. '
                'By default only xyz files are written. '
                'Can only be specified if "ocean_basin_points_filename" is not specified '
                '(ie, ocean basin points must be on a uniform lon/lat grid).')
    
    parser.add_argument('-i', '--ocean_basin_grid_spacing', type=float,
            help='The grid spacing (in degrees) of ocean basin points in lon/lat space. '
                'The grid point latitudes/longitudes are offset by half the grid spacing '
                '(eg, for a 1 degree spacing the latitudes are -89.5, -88.5, ..., 89.5). '
                'Can only be specified if "ocean_basin_points_filename" is not specified. '
                'Defaults to {0} degrees.'.format(
                        DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES))
    
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
                'Can only be specified if "ocean_basin_grid_spacing" and "output_grd_files" are not specified.')
    
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
        not args.output_overriding_plate_with_time and
        not args.output_mean_distance and
        not args.output_std_dev_distance):
        args.output_distance_with_time = True
    
    if (args.output_all_proximity_distances and
        not args.non_topological_proximity_features):
        raise argparse.ArgumentTypeError(
            "Can only specify 'output_all_proximity_distances' if 'non_topological_proximity_features' is also specified.")
    
    # Make sure user provides *topological* features and specifies subduction zone when
    # they want to output overriding plate IDs.
    if args.output_overriding_plate_with_time:
        if args.non_topological_proximity_features:
            raise argparse.ArgumentTypeError(
                "'output_overriding_plate' and 'non_topological_proximity_features' cannot both be specified.")
        if (not args.proximity_feature_types or
            len(args.proximity_feature_types) != 1 or
            args.proximity_feature_types[0] != pygplates.FeatureType.gpml_subduction_zone.get_name()):
            raise argparse.ArgumentTypeError(
                "'proximity_feature_types' must specify only 'SubductionZone' when also specifying 'output_overriding_plate'.")
    
    if args.ocean_basin_points_filename is not None:
        if args.ocean_basin_grid_spacing is not None:
            raise argparse.ArgumentTypeError(
                "'ocean_basin_grid_spacing' and 'ocean_basin_points_filename' cannot both be specified.")
        if args.output_grd_files is not None:
            raise argparse.ArgumentTypeError(
                "'output_grd_files' and 'ocean_basin_points_filename' cannot both be specified.")
    
    # Get the input points.
    if args.ocean_basin_points_filename is not None:
        input_points = read_input_points(args.ocean_basin_points_filename)
    else:
        if args.ocean_basin_grid_spacing is None:
            args.ocean_basin_grid_spacing = DEFAULT_GRID_INPUT_POINTS_GRID_SPACING_DEGREES
        input_points, num_grid_longitudes, num_grid_latitudes = generate_input_points_grid(args.ocean_basin_grid_spacing)
    
    if input_points is None:
        sys.exit(1)
    
    proximity_data = proximity(
            input_points,
            args.rotation_filenames,
            args.proximity_filenames,
            not args.non_topological_proximity_features, # proximity_features_are_topological
            args.proximity_feature_types,
            args.topological_reconstruction_filenames,
            args.max_topological_reconstruction_time,
            args.age_grid_filename,
            args.age_grid_paleo_time,
            args.time_increment,
            args.output_distance_with_time,
            args.output_overriding_plate_with_time,
            args.output_mean_distance,
            args.output_std_dev_distance,
            args.output_all_proximity_distances,
            args.num_cpus,
            args.anchor_plate_id,
            args.max_distance)

    if proximity_data is None:
        sys.exit(1)
    
    write_proximity_data(
            proximity_data,
            args.output_filename_prefix,
            args.output_filename_extension,
            args.time_increment,
            args.output_distance_with_time,
            args.output_overriding_plate_with_time,
            args.output_mean_distance,
            args.output_std_dev_distance,
            (args.ocean_basin_grid_spacing, num_grid_longitudes, num_grid_latitudes) if args.output_grd_files else None)
    
    sys.exit(0)
