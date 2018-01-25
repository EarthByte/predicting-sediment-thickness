# Combine all files into one array

# paste sedthick_1deg.xyz sed_rate_cm_per_ka_1deg.z\
paste average_sedimentation_rate_from_decomp_sed_thick_cm_per_ka_70_80.xyz \
sedthick_1deg.z  \
agegrid_areps_1deg_masked.z  \
PCB_mean_dist_grid_70_80_1deg.z  \
river_mouth_proximity_mean_distance_Amazon_1deg.z  \
river_mouth_proximity_mean_distance_Congo_1deg.z  \
river_mouth_proximity_mean_distance_Ganges-Bramaputra_1deg.z  \
river_mouth_proximity_mean_distance_Godavari-Krishna_1deg.z  \
river_mouth_proximity_mean_distance_Indus_1deg.z  \
river_mouth_proximity_mean_distance_Magdalena_1deg.z \
river_mouth_proximity_mean_distance_Mahanadi_1deg.z \
river_mouth_proximity_mean_distance_Mississippi_1deg.z \
river_mouth_proximity_mean_distance_Narmada_1deg.z \
river_mouth_proximity_mean_distance_Niger_1deg.z \
river_mouth_proximity_mean_distance_Orinoco_1deg.z \
river_mouth_proximity_mean_distance_Paleo-Congo_1deg.z \
river_mouth_proximity_mean_distance_Parana_1deg.z  \
| awk '{if ($3 != "NaN") print}' > alldata.xyz

