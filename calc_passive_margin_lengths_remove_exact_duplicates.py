
import os
import pygplates
import sys

rotation_file = 'D:/Users/john/Downloads/gplates/data/PlateModels/Muller_etal_AREPS_Supplement/Global_EarthByte_230-0Ma_GK07_AREPS.rot'
proximity_features_file = 'Global_EarthByte_GeeK07_COBLineSegments_2016_v3.gpmlz'


#
# Private function.
#
# Removes duplicate reconstructed feature geometries from the list 'reconstructed_feature_geometries'.
# This can reduce the number of reconstructed feature geometries in the list.
# 
def _remove_duplicate_reconstructed_feature_geometries(reconstructed_feature_geometries):
    if len(reconstructed_feature_geometries) == 1:
        return
    
    #
    # Do N^2 search over pairs of reconstructed feature geometries to test for duplicates.
    #
    
    # Iterate over all reconstructed feature geometries except last reconstructed feature geometry.
    # Using len() since some reconstructed feature geometries are removed during iteration.
    rfg1_index = 0
    while rfg1_index < len(reconstructed_feature_geometries) - 1:
        rfg1 = reconstructed_feature_geometries[rfg1_index]
        rfg1_geom = rfg1.get_reconstructed_geometry()
        
        # Iterate over the remaining reconstructed feature geometries (after 'rfg1').
        # Using len() since some reconstructed feature geometries are removed during iteration.
        rfg2_index = rfg1_index + 1
        while rfg2_index < len(reconstructed_feature_geometries):
            rfg2 = reconstructed_feature_geometries[rfg2_index]
            rfg2_geom = rfg2.get_reconstructed_geometry()
            
            # Compare the geometries of rfg1 and rfg2.
            # Test for duplicate geometries.
            if rfg1_geom == rfg2_geom:
                del reconstructed_feature_geometries[rfg2_index]
                rfg2_index -= 1
            
            rfg2_index += 1
        
        rfg1_index += 1


rotation_model = pygplates.RotationModel(rotation_file)

proximity_features = list(pygplates.FeatureCollection(proximity_features_file))

print('{:<10} {:<10}'.format('Age(Ma)', 'Length(Kms)'))

# Time range 0-230Ma in 1My intervals.
for time in range(0, 231):
    
    length_in_kms = 0
    
    # Reconstruct the proximity features that exist at the current 'time'.
    proximity_reconstructed_feature_geometries = []
    pygplates.reconstruct(proximity_features, rotation_model, proximity_reconstructed_feature_geometries, time)
    
    # Remove any 'exact' duplicate reconstructed geometries.
    # This does not remove partially overlapping geometries, etc.
    _remove_duplicate_reconstructed_feature_geometries(proximity_reconstructed_feature_geometries)
    
    # Iterate over reconstructed geometries and add up their segment lengths.
    for proximity_reconstructed_feature_geometry in proximity_reconstructed_feature_geometries:
        proximity_reconstructed_geometry = proximity_reconstructed_feature_geometry.get_reconstructed_geometry()
        
        try:
            length_in_kms += proximity_reconstructed_geometry.get_arc_length() * pygplates.Earth.mean_radius_in_kms
        except AttributeError:
            # Only polylines and polygons have the 'get_arc_length()' method.
            # Skip points and multipoints...
            continue
    
    print('{:<10.0f} {:<10.2f}'.format(time, length_in_kms))
