from __future__ import print_function
import os
import pygplates
import sys

rotation_file = 'D:/Users/john/Downloads/gplates/data/PlateModels/Muller_etal_AREPS_Supplement/Global_EarthByte_230-0Ma_GK07_AREPS.rot'
proximity_features_file = 'Global_EarthByte_GeeK07_COBLineSegments_2016_v3.gpmlz'


rotation_model = pygplates.RotationModel(rotation_file)
proximity_features = list(pygplates.FeatureCollection(proximity_features_file))

print('{:<10} {:<10}'.format('Age(Ma)', 'Length(Kms)'))

# Time range 0-230Ma in 1My intervals.
for time in range(0, 231):
    
    length_in_kms = 0
    
    # Reconstruct the proximity features that exist at the current 'time'.
    proximity_reconstructed_feature_geometries = []
    pygplates.reconstruct(proximity_features, rotation_model, proximity_reconstructed_feature_geometries, time)
    
    for proximity_reconstructed_feature_geometry in proximity_reconstructed_feature_geometries:
        proximity_reconstructed_geometry = proximity_reconstructed_feature_geometry.get_reconstructed_geometry()
        
        try:
            length_in_kms += proximity_reconstructed_geometry.get_arc_length() * pygplates.Earth.mean_radius_in_kms
        except AttributeError:
            # Only polylines and polygons have the 'get_arc_length()' method.
            # Skip points and multipoints...
            continue
    
    print('{:<10.0f} {:<10.2f}'.format(time, length_in_kms))
