# -*- coding: utf-8 -*-

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


#######################################################################
# Efficient point-to-geometry distance queries when there are many    #
# relatively uniformly spaced points to be tested against geometries. #
#######################################################################
#
#
# For example, to find the closest geometry to each point in a sequence of points:
#
#
#    import distance_points_to_geometries
#    import math
#    
#    # A list of 'pygplates.PointOnSphere' points.
#    points = [...]
#
#    # Some geometry features (eg, coastlines).
#    geometry_feature_collection = pygplates.FeatureCollection('geometries.gpml')
#
#    # Look for features within 90 degrees of each point.
#    distance_threshold_radians = math.pi / 2
#
#    # Extract the geometries from the features.
#    geometries = []
#    geometry_features = []
#    for geometry_feature in geometry_feature_collection:
#        geometries.append(geometry_feature.get_geometry())
#        geometry_features.append(geometry_feature)
#
#    #
#    # Find the closest geometry (feature) to each point (within threshold distance).
#    #
#    geometry_features_closest_to_points = distance_points_to_geometries.distances(
#            points,
#            geometries,
#            geometry_features,
#            distance_threshold_radians = distance_threshold_radians)
#
#    # Print name of closest feature to each point (if any).
#    for point_index, closest_geometry_feature in enumerate(geometry_features_closest_to_points):
#        if closest_geometry_feature is not None:
#            distance, geometry_feature = closest_geometry_feature
#            print('Closest feature to', points[point_index].to_lat_lon(), 'is', geometry_feature.get_name(),
#                    'with distance', distance * pygplates.Earth.mean_radius_in_kms, 'kms')
#        else:
#            print('No features close to', points[point_index].to_lat_lon())
#
#    #
#    # Find all geometries (features) near each point (within threshold distance).
#    #
#    geometry_features_closest_to_points = distance_points_to_geometries.distances(
#            points,
#            geometries,
#            geometry_features,
#            distance_threshold_radians = distance_threshold_radians,
#            all_geometries=True)
#
#    # Print names of the closest features to each point (if any).
#    for point_index, geometry_feature_list in enumerate(geometry_features_closest_to_points):
#        if geometry_feature_list:
#            print('Closest features to', points[point_index].to_lat_lon(), 'are...')
#            for distance, geometry_feature in geometry_feature_list
#                print('    ', geometry_feature.get_name(), 'with distance',
#                        distance * pygplates.Earth.mean_radius_in_kms, 'kms')
#        else:
#            print('No features close to', points[point_index].to_lat_lon())
#
#####################################################################


from __future__ import print_function
import points_spatial_tree
import math
import pygplates
import sys


def distances(
        points,
        geometries,
        geometry_proxies = None,
        distance_threshold_radians = None,
        return_closest_position = False,
        return_closest_index = False,
        geometries_are_solid = False,
        all_geometries = False,
        subdivision_depth = points_spatial_tree.DEFAULT_SUBDIVISION_DEPTH):
    """
    Efficient point-to-geometry distance queries when there are many relatively uniformly spaced points to be tested against geometries.
    
    points: a sequence of 'pygplates.PointOnSphere'.
    
    geometries: a sequence of 'pygplates.GeometryOnSphere'.
    
    geometry_proxies: Optional sequence of objects associated with 'geometries'.
                     If not specified then the proxies default to the geometries themselves.
                     These can be any object (such as the 'pygplates.Feature' that the geometry came from).
    
    distance_threshold_radians: Optional distance threshold in radians - threshold should be in the range [0,PI] if specified.
    
    return_closest_position: Whether to also return the closest point on each geometry - default is False.
    
    return_closest_index: Whether to also return the index of the closest point (for multi-points) or
                          the index of the closest segment (for polylines and polygons) - default is False.
    
    geometries_are_solid: Whether the interiors of the geometries are solid or not - only applies to polygon geometries - default is False.
    
    all_geometries: Whether to find all geometries near each point (within threshold distance) or just the closest.
                    Defaults to False (only returns closest geometry to each point).
    
    subdivision_depth: The depth of the lat/lon quad tree used to speed up point-to-geometry distance queries.
                       The lat/lon width of a leaf quad tree node is (90 / (2^subdivision_depth)) degrees.
                       Generally the denser the 'points' the larger the depth should be.
                       Setting this value too high causes unnecessary time to be spent generating a deep quad tree.
                       Setting this value too low reduces the culling efficiency of the quad tree.
                       However a value of 4 seems to work quite well for a uniform lat/lon spacing of 'points' of 1 degree and below
                       without the cost of generating a deep quad tree.
                       So most of the time the subdivision depth can be left at its default value.
    
    Returns: A list of geometry proxies associated with 'points'.
             The length of the returned list matches the length of 'points'.
             For each point in 'points', if the point is close to a geometry then that geometry's proxy (and its distance information)
             is stored (otherwise None is stored) at the same index (as the point) in the returned list.
             If 'all_geometries' is False then each item in returned list is a single geometry proxy (and its distance information)
             representing the closest geometry within threshold distance (or a single None).
             If 'all_geometries' is True then each item in returned list is a *list* of geometry proxies (and their distance informations)
             representing all geometries within threshold distance (or a single None).
             Above we mentioned "geometry proxy (and its distance information)". This is a tuple whose size depends on
             the values of 'return_closest_position' and 'return_closest_index' according to...
             
                if return_closest_position and return_closest_index:
                    geometry_proxy_to_point = (distance, closest_position, closest_index, geometry_proxy)
                elif return_closest_position:
                    geometry_proxy_to_point = (distance, closest_position, geometry_proxy)
                elif return_closest_index:
                    geometry_proxy_to_point = (distance, closest_index, geometry_proxy)
                else:
                    geometry_proxy_to_point = (distance, geometry_proxy)
    
    The arguments 'distance_threshold_radians', 'return_closest_position', 'return_closest_index' and 'geometries_are_solid' are
    similar to those in pygplates.GeometryOnSphere.distance()...
    See http://www.gplates.org/docs/pygplates/generated/pygplates.GeometryOnSphere.html#pygplates.GeometryOnSphere.distance
    
    Raises ValueError if the lengths of 'geometries' and 'geometry_proxies' (if specified) do not match.
    """
    
    spatial_tree_of_points = points_spatial_tree.PointsSpatialTree(points, subdivision_depth)
    
    return distances_using_points_spatial_tree(
            points,
            spatial_tree_of_points,
            geometries,
            geometry_proxies,
            distance_threshold_radians,
            return_closest_position,
            return_closest_index,
            geometries_are_solid,
            all_geometries)


def distances_using_points_spatial_tree(
        points,
        spatial_tree_of_points,
        geometries,
        geometry_proxies = None,
        distance_threshold_radians = None,
        return_closest_position = False,
        return_closest_index = False,
        geometries_are_solid = False,
        all_geometries = False):
    """
    Same as 'distances()' except 'spatial_tree_of_points' is a 'points_spatial_tree.PointsSpatialTree' of 'points'.
    
    This is useful when re-using a single 'points_spatial_tree.PointsSpatialTree'.
    For example, when using it both for point-in-polygon queries and minimum distance queries.
    """
    
    # Use the geometries as proxies if no proxies have been specified.
    if geometry_proxies is None:
        geometry_proxies = geometries
    
    if len(geometries) != len(geometries):
        raise ValueError('Number of geometries must match number of geometry proxies.')
    
    # 'geometries_and_proxies' is a list of 2-tuples (geometry, geometry_proxy).
    geometries_and_proxies = [(geometries[index], geometry_proxies[index]) for index in xrange(len(geometries))]
    
    # By default no points are within threshold distance to any geometry.
    # If any are found to be within threshold distance then we'll set the proxy of the closest geometry.
    geometry_proxies_closest_to_points = [None] * len(points)
    
    for root_node in spatial_tree_of_points.get_root_nodes():
        _visit_spatial_tree_node(
                root_node,
                points,
                geometries_and_proxies,
                geometry_proxies_closest_to_points,
                distance_threshold_radians,
                return_closest_position,
                return_closest_index,
                geometries_are_solid,
                all_geometries)
    
    return geometry_proxies_closest_to_points


##################
# Implementation #
##################


def _visit_spatial_tree_node(
        node,
        points,
        parent_closest_geometries_and_proxies,
        geometry_proxies_closest_to_points,
        distance_threshold_radians,
        return_closest_position,
        return_closest_index,
        geometries_are_solid,
        all_geometries):
    
    node_bounding_circle_centre, node_bounding_circle_radius = node.get_bounding_circle()
        
    # If there is a distance threshold then ignore the geometries that are further
    # from node bounding circle than the distance threshold.
    # Note that the distance to the bounding circle is the distance to its *centre* minus its radius.
    # In other words, the distance to its *centre* is the distance to the circle plus its radius.
    distance_threshold_to_node_centre_radians = distance_threshold_radians
    if distance_threshold_to_node_centre_radians is not None:
        distance_threshold_to_node_centre_radians += node_bounding_circle_radius
        # If threshold exceeds maximum possible distance then we don't need a threshold.
        if distance_threshold_to_node_centre_radians >= math.pi:
            distance_threshold_to_node_centre_radians = None
    
    # See if the current quad tree node's bounding circle is close to any geometries.
    closest_geometries_and_proxies = []
    
    # Some extra parameters if we're only interested in the *closest* geometry to each point.
    if not all_geometries:
        distance_node_centre_to_geometries = []
        min_distance_node_centre_to_geometries = None
        max_distance_node_centre_to_geometries = None
    for geometry, geometry_proxy in parent_closest_geometries_and_proxies:
            
        distance_node_centre_to_geometry = pygplates.GeometryOnSphere.distance(
                node_bounding_circle_centre,
                geometry,
                distance_threshold_to_node_centre_radians,
                geometry2_is_solid = geometries_are_solid)
        if distance_node_centre_to_geometry is None:
            # Ignore the current geometry.
            # This can only happen if the distance threshold is not None.
            continue
        
        closest_geometries_and_proxies.append((geometry, geometry_proxy))
        
        # Calculate extra parameters if we're only interested in the *closest* geometry to each point.
        if not all_geometries:
            distance_node_centre_to_geometries.append(
                    (distance_node_centre_to_geometry, len(closest_geometries_and_proxies) - 1))
            
            if (min_distance_node_centre_to_geometries is None or
                distance_node_centre_to_geometry < min_distance_node_centre_to_geometries):
                
                min_distance_node_centre_to_geometries = distance_node_centre_to_geometry
            
            if (max_distance_node_centre_to_geometries is None or
                distance_node_centre_to_geometry > max_distance_node_centre_to_geometries):
                
                max_distance_node_centre_to_geometries = distance_node_centre_to_geometry
    
    # If quad tree node is further than threshold distance to all geometries then nothing to do since
    # all points are marked as not within the distance threshold.
    if not closest_geometries_and_proxies:
        return
    
    # If we're only interested in the *closest* geometry to each point then we can exclude
    # geometries that cannot possibly be the closest to any point in the current node.
    if not all_geometries:
        # If the difference between the min and max distances exceeds the bounding circle diameter then
        # we can remove those geometries whose minimum distance to bounding circle exceeds the maximum
        # distance to the closest geometry - this means the distance from the closest geometry to the
        # point (in the current node) *furthest* from it is still less than the distance from another geometry
        # to the point (in the current node) *closest* to that geometry - hence all points (in the current node)
        # are closer to the closest geometry.
        node_bounding_circle_diameter = 2 * node_bounding_circle_radius
        if max_distance_node_centre_to_geometries - min_distance_node_centre_to_geometries > node_bounding_circle_diameter:
            # Sort from smallest to largest distance to make geometry removal easier.
            distance_node_centre_to_geometries.sort()
            
            # Essentially remove any geometries whose distance compared to the closest geometry exceeds
            # the diameter of the current node's bounding circle.
            new_closest_geometries_and_proxies = []
            for index in xrange(0, len(distance_node_centre_to_geometries)):
                distance_node_centre_to_geometry, closest_geometries_and_proxies_index = distance_node_centre_to_geometries[index]
                if distance_node_centre_to_geometry - min_distance_node_centre_to_geometries > node_bounding_circle_diameter:
                    # All remaining geometries are essentially removed since the distance list is sorted.
                    break
                new_closest_geometries_and_proxies.append(
                        closest_geometries_and_proxies[closest_geometries_and_proxies_index])
            
            closest_geometries_and_proxies = new_closest_geometries_and_proxies
    
    # Visit child nodes (if internal node) or test each point (if leaf node).
    if node.is_internal_node():
        for child_node in node.get_child_nodes():
            _visit_spatial_tree_node(
                    child_node,
                    points,
                    closest_geometries_and_proxies,
                    geometry_proxies_closest_to_points,
                    distance_threshold_radians,
                    return_closest_position,
                    return_closest_index,
                    geometries_are_solid,
                    all_geometries)
    else:
        for point_index in node.get_point_indices():
            point = points[point_index]
            
            # Start out with the distance threshold (which might be None).
            # If only looking for closest geometry then we'll reduce this to the closest geometry so far as we go.
            distance_threshold_to_point = distance_threshold_radians
            
            # If only looking for closest geometry then keep track of it.
            if not all_geometries:
                closest_geometry_proxy_to_point = None
            
            for geometry, geometry_proxy in closest_geometries_and_proxies:
                
                point_to_geometry_distance_info = pygplates.GeometryOnSphere.distance(
                        point,
                        geometry,
                        distance_threshold_to_point,
                        return_closest_position,
                        return_closest_index,
                        # Whether to treat 'geometry' as solid or not (if it's a polygon)...
                        geometry2_is_solid = geometries_are_solid)
                
                # If point is close to a geometry (or closer than previous closest geometry).
                if point_to_geometry_distance_info is not None:
                    # Unpack the distance info.
                    if return_closest_position:
                        if return_closest_index:
                            distance, _, closest_position, _, closest_index = point_to_geometry_distance_info
                            geometry_proxy_to_point = (distance, closest_position, closest_index, geometry_proxy)
                        else:
                            distance, _, closest_position = point_to_geometry_distance_info
                            geometry_proxy_to_point = (distance, closest_position, geometry_proxy)
                    elif return_closest_index:
                        distance, _, closest_index = point_to_geometry_distance_info
                        geometry_proxy_to_point = (distance, closest_index, geometry_proxy)
                    else:
                        distance = point_to_geometry_distance_info
                        geometry_proxy_to_point = (distance, geometry_proxy)
                    
                    if all_geometries:
                        # Each point has a *list* of geometry proxies (or None).
                        # Create list if first geometry proxy encountered for current point.
                        if geometry_proxies_closest_to_points[point_index] is None:
                            geometry_proxies_closest_to_points[point_index] = []
                        geometry_proxies_closest_to_points[point_index].append(geometry_proxy_to_point)
                    else:
                        # Only looking for closest geometry, so reduce distance threshold to closest so far.
                        distance_threshold_to_point = distance
                        closest_geometry_proxy_to_point = geometry_proxy_to_point
            
            # If only looking for closest geometry then only need to store closest geometry (not a list of geometries).
            if not all_geometries:
                if closest_geometry_proxy_to_point is not None:
                    geometry_proxies_closest_to_points[point_index] = closest_geometry_proxy_to_point


#if __name__ == '__main__':
#    
#    #
#    # Some testing/example code.
#    #
#    
#    import time
#    
#    
#    print('Loading coastline polygons and rotation model...')
#    coastline_features = pygplates.FeatureCollection('../../../sample_data/2.0/SampleData/FeatureCollections/Coastlines/Matthews_etal_GPC_2016_Coastlines.gpmlz')
#    rotation_model = pygplates.RotationModel('../../../sample_data/2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot')
#    
#    print('Reconstructing coastline polygons...')
#    reconstruction_time = 200
#    coastline_reconstructed_feature_geometries = []
#    pygplates.reconstruct(coastline_features, rotation_model, coastline_reconstructed_feature_geometries, reconstruction_time)
#    
#    geometries = []
#    geometry_features = []
#    for reconstructed_feature_geometry in coastline_reconstructed_feature_geometries:
#        geometries.append(reconstructed_feature_geometry.get_reconstructed_geometry())
#        geometry_features.append(reconstructed_feature_geometry.get_feature())
#    
#    # Create uniform lat/lon distribution of points.
#    print('Creating lat/lon grid of points...')
#    num_latitudes = 180
#    num_longitudes = 360
#    lat_grid_spacing_degrees = 180.0 / num_latitudes
#    lon_grid_spacing_degrees = 360.0 / num_longitudes
#    
#    points = []
#    for lat_index in xrange(num_latitudes):
#        # The 0.5 puts the point in the centre of the grid pixel.
#        # This also avoids sampling right on the poles.
#        lat = -90 + (lat_index + 0.5) * lat_grid_spacing_degrees
#        
#        for lon_index in xrange(num_longitudes):
#            # The 0.5 puts the point in the centre of the grid pixel.
#            # This also avoids sampling right on the dateline where there might be
#            # age grid or static polygon artifacts.
#            lon = -180 + (lon_index + 0.5) * lon_grid_spacing_degrees
#            
#            point = pygplates.PointOnSphere(lat, lon)
#            points.append(point)
#    
#    print('Finding geometries closest to points...')
#    time_begin = time.clock()
#    
#    distance_threshold_radians = 400 / pygplates.Earth.mean_radius_in_kms
#    
#    if True:
#        #
#        # The fast way (about 3 seconds).
#        #
#        geometry_features_closest_to_points = distances(
#                points,
#                geometries,
#                geometry_features,
#                distance_threshold_radians = distance_threshold_radians,
#                all_geometries=True)
#    else:
#        #
#        # The slow way (about 270 seconds).
#        #
#        # Similar to 'distances()' except without using a quad tree.
#        #
#        geometries_and_features = [(geometries[index], geometry_features[index]) for index in xrange(len(geometries))]
#        geometry_features_closest_to_points = [None] * len(points)
#        for point_index, point in enumerate(points):
#            geometry_features_closest_to_points[point_index] = []
#            for geometry, geometry_feature in geometries_and_features:
#                dist = pygplates.GeometryOnSphere.distance(point, geometry, distance_threshold_radians)
#                if dist is not None:
#                    geometry_features_closest_to_points[point_index].append((dist, geometry_feature))
#    
#    time_end = time.clock()
#    print('  {0} seconds'.format(time_end - time_begin))
#    
#    print('Associate each point with zero or more closest geometries...')
#    
#    # Group points with each geometry so can create one multi-point per geometry.
#    geometry_feature_to_points_mapping = {}
#    for point_index, geometry_feature_list in enumerate(geometry_features_closest_to_points):
#        if geometry_feature_list:
#            for distance, geometry_feature in geometry_feature_list:
#                points_near_geometry, distances_near_geometry = geometry_feature_to_points_mapping.setdefault(geometry_feature, ([], []))
#                points_near_geometry.append(points[point_index])
#                distances_near_geometry.append(distance * pygplates.Earth.mean_radius_in_kms)
#    
#    # Create multi-point features.
#    multi_point_features = []
#    for geometry_feature, (points_near_geometry, distances_near_geometry) in geometry_feature_to_points_mapping.iteritems():
#        multi_point_feature = pygplates.Feature()
#        multi_point_feature.set_geometry((
#                pygplates.MultiPointOnSphere(points_near_geometry),
#                {pygplates.ScalarType.create_gpml('Distance') : distances_near_geometry}))
#        
#        begin_time, end_time = geometry_feature.get_valid_time()
#        multi_point_feature.set_valid_time(begin_time, end_time)
#        
#        multi_point_feature.set_reconstruction_plate_id(
#                geometry_feature.get_reconstruction_plate_id())
#        
#        multi_point_features.append(multi_point_feature)
#    
#    print('Writing points feature collection...')
#    pygplates.FeatureCollection(multi_point_features).write('multi_point_features.gpml')
