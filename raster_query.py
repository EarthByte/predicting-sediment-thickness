
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


################################################################################################
# Efficiently query values in rasters/grids at point locations including:                      #
# - Query at arbitrary points.                                                                 #
# - Query at tessellated points along reconstructed geometries or resolved topologies.         #
#   Note: Tessellation is only needed for polylines and polygons (not points and multipoints). #
################################################################################################


from __future__ import print_function
import argparse
from call_system_command import call_system_command
import math
import points_spatial_tree
import proximity_query
import pygplates
import sys


def query_raster_at_points(
        raster_filename,
        query_points,
        search_radius_radians,
        smoothing_radius_radians = None):
    """
    Query the raster at query points.
    
    raster_filename: The filename of the raster/grid file.
    
    query_points: The query points as a sequence of pygplates.PointOnSphere.
    
    search_radius_radians: The distance (in radians) from query point to search for non-NaN raster grid values.
                           If none are found for a query point then it will have a query scalar value of NaN.
    
    smoothing_radius_radians: Determines which non-NaN raster grid values (if any) to use in a weighted average for a query point.
                              All points within a radius of 'min_dist + smoothing_radius_radians' of the query point are included
                              (where 'min_dist' is the distance to the closest non-NaN raster grid location).
                              Note that 'smoothing_radius_radians' should be less than 'search_radius_radians'.
                              Default value is None which results in only the closest non-NaN raster grid value being chosen.
    
    Returns a list raster scalars associated with the query points (one scalar per query point).
    Each scalar value is either:
    - sampled from raster (at query point location), or
    - the nearest raster value(s) if query point is in masked (NaN) region of raster (and within search radius), or
    - NaN if there are no non-NaN grid values within search radius.
    """
    
    # Get a list of points and a list of scalar for each pixel node location in raster (excluding no-data nodes).
    raster_grid_points, raster_grid_scalars = get_raster_grid_positions_and_scalars(raster_filename)
    
    # Reuse the spatial tree of raster grid points since it's expensive to generate.
    raster_grid_points_spatial_tree = points_spatial_tree.PointsSpatialTree(raster_grid_points)
    
    # Find the closest raster grid points near each query point (within threshold distance).
    #
    # Note that this only finds the closest raster grid point to each query point.
    # If smoothing is enabled then another search will be done with a smoothing radius around closest point (per-query-point).
    raster_grid_points_closest_to_query_points = proximity_query.find_closest_points_to_geometries_using_points_spatial_tree(
            query_points,
            raster_grid_points,
            raster_grid_points_spatial_tree,
            # Indices into the raster grid points and scalars...
            range(len(raster_grid_points)),
            distance_threshold_radians = search_radius_radians)
    
    # List of associated raster values at query point locations.
    query_scalars = []
    
    # Iterate over the query points to get the closest raster grid point (to each query point).
    for query_point_index, closest_raster_grid_point_distance_and_index in enumerate(raster_grid_points_closest_to_query_points):
        if closest_raster_grid_point_distance_and_index is not None:
            # Get the closest raster grid point to the current query point.
            _, closest_raster_grid_point_index = closest_raster_grid_point_distance_and_index
            
            # We'll need to collect more than just the closest raster grid value if smoothing is enabled.
            if smoothing_radius_radians is not None:
                # Find the raster grid points within smoothing radius of closest raster grid point to the current query point.
                raster_grid_point_scalars_to_weight = proximity_query.find_closest_points_to_geometry_using_points_spatial_tree(
                        raster_grid_points[closest_raster_grid_point_index],
                        raster_grid_points,
                        raster_grid_points_spatial_tree,
                        raster_grid_scalars,
                        distance_threshold_radians = smoothing_radius_radians,
                        all_points=True)
                
                # Calculate distance-weighted average of closest raster scalar values.
                sum_weights = 0.0
                sum_weighted_scalars = 0.0
                for distance, scalar in raster_grid_point_scalars_to_weight:
                    # Weight range is [0.1, 1.0] for distance range [smoothing_radius_radians, 0.0].
                    weight = ((smoothing_radius_radians * smoothing_radius_radians) /
                            (smoothing_radius_radians * smoothing_radius_radians + 9 * distance * distance))
                    sum_weights += weight
                    sum_weighted_scalars += weight * scalar
                
                query_point_scalar = sum_weighted_scalars / sum_weights
            else:
                # No smoothing so just choose the scalar value of the closest raster grid point.
                query_point_scalar = raster_grid_scalars[closest_raster_grid_point_index]
        
        else:
            query_point_scalar = float('nan')
        
        query_scalars.append(query_point_scalar)
    
    return query_scalars


def query_raster_with_reconstructed_geometries(
        raster_filename,
        rotation_features_or_model,
        time,
        query_features,
        tessellation_threshold_radians,
        search_radius_radians,
        smoothing_radius_radians = None,
        query_feature_types = None,
        anchor_plate_id = 0):
    """
    Reconstructs and tessellates regular features, and queries raster at the reconstructed tessellated positions.
    
    raster_filename: The filename of the raster/grid file.
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    time: The reconstruction time.
    
    query_features: Regular feature collection(s), or list of features, or filename(s) or any combination of those.
    
    tessellation_threshold_radians: Threshold sampling distance along polylines/polygons (in radians).
    
    search_radius_radians: The distance (in radians) from query point to search for non-NaN raster grid values.
                           If none are found for a query point then it will have a query scalar value of NaN.
    
    smoothing_radius_radians: Determines which non-NaN raster grid values (if any) to use in a weighted average for a query point.
                              All points within a radius of 'min_dist + smoothing_radius_radians' of the query point are included
                              (where 'min_dist' is the distance to the closest non-NaN raster grid location).
                              Note that 'smoothing_radius_radians' should be less than 'search_radius_radians'.
                              Default value is None which results in only the closest non-NaN raster grid value being chosen.
    
    query_feature_types: Optional sequence of feature types to filter/accept (defaults to all feature types).
    
    Returns a list with a tuple for each query point containing the following parameters:
    - query point longitude
    - query point latitude
    - query scalar value, is either:
      + sampled from raster (at query point location), or
      + the nearest raster value(s) if query point is in masked (NaN) region of raster (and within search radius), or
      + NaN if there are no non-NaN grid values within search radius.
    - length of tessellated arc segment (in radians) that query point is on (if a tessellated point on polyline/polygon), or
      NaN if query point is from a point or multipoint geometry.
    
    Only polylines and polygons are tessellated. Points and multipoints just return their points.
    Each tessellated point is the midpoint of a segment of the tessellated polyline/polygon.
    The tessellated weight is the arc length (in radians) of the associated segments.
    """
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn query data into a list of features (if not already).
    query_features = pygplates.FeaturesFunctionArgument(query_features).get_features()
    
    if query_feature_types:
        # Remove those features not matching the allowed feature types.
        # Create a new list containing only features matching the allowed feature types.
        query_features = [feature for feature in query_features
                if feature.get_feature_type() in query_feature_types]
    
    # Reconstruct features that exist at the current 'time'.
    query_reconstructed_feature_geometries = []
    pygplates.reconstruct(query_features, rotation_model, query_reconstructed_feature_geometries, time, anchor_plate_id)
    
    query_reconstructed_geometries = []
    for query_reconstructed_feature_geometry in query_reconstructed_feature_geometries:
        query_reconstructed_geometries.append(query_reconstructed_feature_geometry.get_reconstructed_geometry())
    
    # Tessellate query geometries to get query points and weights.
    query_points, query_point_weights = tessellate_geometries(query_reconstructed_geometries, tessellation_threshold_radians)
    
    # Query the raster at the query points.
    query_point_scalars = query_raster_at_points(
            raster_filename,
            query_points,
            search_radius_radians,
            smoothing_radius_radians)
    
    query_data = []
    for query_point_index in xrange(len(query_points)):
        query_point_lat, query_point_lon = query_points[query_point_index].to_lat_lon()
        query_point_scalar = query_point_scalars[query_point_index]
        query_point_weight = query_point_weights[query_point_index]
        
        # The data will be output in GMT format (ie, lon first, then lat, etc).
        query_data.append((
                query_point_lon,
                query_point_lat,
                query_point_scalar,
                query_point_weight))
    
    return query_data


def query_raster_with_resolved_topologies(
        raster_filename,
        rotation_features_or_model,
        time,
        query_features,
        tessellation_threshold_radians,
        search_radius_radians,
        smoothing_radius_radians = None,
        query_feature_types = None,
        anchor_plate_id = 0):
    """
    Resolves and tessellates topological features, and queries raster at the resolved tessellated positions.
    
    raster_filename: The filename of the raster/grid file.
    
    rotation_features_or_model: Rotation model or feature collection(s), or list of features, or filename(s).
    
    time: The reconstruction time.
    
    query_features: Topological feature collection(s), or list of features, or filename(s) or any combination of those.
    
    tessellation_threshold_radians: Threshold sampling distance along resolved topological section polylines (in radians).
    
    search_radius_radians: The distance (in radians) from query point to search for non-NaN raster grid values.
                           If none are found for a query point then it will have a query scalar value of NaN.
    
    smoothing_radius_radians: Determines which non-NaN raster grid values (if any) to use in a weighted average for a query point.
                              All points within a radius of 'min_dist + smoothing_radius_radians' of the query point are included
                              (where 'min_dist' is the distance to the closest non-NaN raster grid location).
                              Note that 'smoothing_radius_radians' should be less than 'search_radius_radians'.
                              Default value is None which results in only the closest non-NaN raster grid value being chosen.
    
    query_feature_types: Optional sequence of feature types to filter/accept (defaults to all feature types).
                         Note: The feature types apply to the topological *sections* (eg, subduction zone).
    
    Returns a list with a tuple for each query point containing the following parameters:
    - query point longitude
    - query point latitude
    - query scalar value, is either:
      + sampled from raster (at query point location), or
      + the nearest raster value(s) if query point is in masked (NaN) region of raster (and within search radius), or
      + NaN if there are no non-NaN grid values within search radius.
    - length of tessellated polyline arc segment (in radians) that query point is on.
    """
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn query data into a list of features (if not already).
    query_features = pygplates.FeaturesFunctionArgument(query_features).get_features()
    
    # Note that, for *topological* features, we cannot remove those not matching the allowed feature types
    # because we need to resolve *all* topologies and then filter out the shared topology sections by feature type.
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    query_resolved_topologies = []
    query_shared_boundary_sections = []
    pygplates.resolve_topologies(query_features, rotation_model, query_resolved_topologies, time, query_shared_boundary_sections, anchor_plate_id)
    
    # Iterate over the shared boundary sections of all resolved topologies.
    query_reconstructed_geometries = []
    for query_shared_boundary_section in query_shared_boundary_sections:
        # Skip sections that are not included in the list of boundary feature types (if any).
        query_feature = query_shared_boundary_section.get_feature()
        if (query_feature_types and
            query_feature.get_feature_type() not in query_feature_types):
            continue
        
        # Iterate over the shared sub-segments of the current boundary line.
        # These are the parts of the boundary line that actually contribute to topological boundaries.
        for query_shared_sub_segment in query_shared_boundary_section.get_shared_sub_segments():
            query_reconstructed_geometries.append(query_shared_sub_segment.get_resolved_geometry())
    
    # Tessellate query geometries to get query points and weights.
    query_points, query_point_weights = tessellate_geometries(query_reconstructed_geometries, tessellation_threshold_radians)
    
    # Query the raster at the query points.
    query_point_scalars = query_raster_at_points(
            raster_filename,
            query_points,
            search_radius_radians,
            smoothing_radius_radians)
    
    query_data = []
    for query_point_index in xrange(len(query_points)):
        query_point_lat, query_point_lon = query_points[query_point_index].to_lat_lon()
        query_point_scalar = query_point_scalars[query_point_index]
        query_point_weight = query_point_weights[query_point_index]
        
        # The data will be output in GMT format (ie, lon first, then lat, etc).
        query_data.append((
                query_point_lon,
                query_point_lat,
                query_point_scalar,
                query_point_weight))
    
    return query_data


def tessellate_geometries(
        query_geometries,
        tessellation_threshold_radians):
    """
    Tessellate query geometries to get query points and weights.
    
    query_geometries: A sequence of pygplates.GeometryOnSphere (ie, point/multipoint/polyline/polygon).
    
    tessellation_threshold_radians: The maximum spacing of tessellated points along polylines/polygons.
    
    Returns: A 2-tuple of lists.
             The first list contains tessellated pygplates.PointOnSphere's.
             The second list contains the weight as the arc length (in radians) of a tessellated
             segment of polyline/polygon, or NaN for a point (ie, a point geometry or points in a multipoint).
    
    Only polylines and polygons are tessellated. Points and multipoints just return their points.
    Each tessellated point is the midpoint of a segment of the tessellated polyline/polygon.
    The tessellated weight is the arc length (in radians) of the associated segments.
    """
    
    query_points = []
    query_point_weights = []
    
    # Ensure each geometry is tessellated to within the threshold sampling distance.
    for query_geometry in query_geometries:
        try:
            tessellated_query_geometry = query_geometry.to_tessellated(tessellation_threshold_radians)
            
            # Iterate over the great circle arcs of the tessellated polyline/polygon to get the arc midpoints and lengths.
            # There is an arc between each adjacent pair of points in the polyline/polygon.
            for arc in tessellated_query_geometry.get_segments():
                if not arc.is_zero_length():
                    query_points.append(arc.get_arc_point(0.5))
                    query_point_weights.append(arc.get_arc_length())
            
        except AttributeError:
            # Geometry is a point or multipoint - just leave it as is (not need to tessellate).
            query_geometry_points = query_geometry.get_points()
            query_points.extend(query_geometry_points)
            query_point_weights.extend([float('nan')] * len(query_geometry_points))
    
    return query_points, query_point_weights


def get_raster_grid_positions_and_scalars(raster_filename, nodata_value=None):
    """
    Returns a 2-tuple of lists.
    The first is a list of pygplates.PointOnSphere at all pixel node locations in raster.
    The second is a list of raster scalar values at those pixel node locations.
    Both list have the same length.
    
    Note that this excludes points that have the no-data value.
    If 'nodata_value' is not specified then it defaults to NaN.
    """

    # The command-line strings to execute GMT 'grd2xyz'.
    grd2xyz_command_line = ["gmt", "grd2xyz", "-s"]
    if nodata_value is not None:
        grd2xyz_command_line.append("-di{0}".format(nodata_value))
    grd2xyz_command_line.append(raster_filename)
    
    stdout_data = call_system_command(grd2xyz_command_line, return_stdout=True)
    
    # Create a list of points and a list of scalars.
    raster_grid_points = []
    raster_grid_scalars = []
    
    # Read lon, lat and scalar values from the output of 'grd2xyz'.
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
            
            # The scalar got appended to the last column by 'grd2xyz'.
            # Note: We should not have NaN values because the '-s' option to 'grd2xyz' removes them.
            scalar = float(line_data[-1])
            
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point lon, lat and scalar values.'.format(line), file=sys.stderr)
            continue
        
        raster_grid_points.append(pygplates.PointOnSphere(lat, lon))
        raster_grid_scalars.append(scalar)
    
    return raster_grid_points, raster_grid_scalars
