#!/bin/zsh

# mean distance to passive cont margin grid

# i=PCB_mean_dist_grid_70_80.grd
i=mean_distance_1d_0.grd
j=`basename $i .grd`

if [ ! -s ${j}_1deg.z ]; then
	echo "Converting $i to $j.z"
	grdsample -I1 -fg -R-180/180/-70/80 $i -G${j}_1deg.grd
	grd2xyz ${j}_1deg.grd | cut -f3 >${j}_1deg.z
fi

exit

# convert sed thickness grid

# decompacted sed rate grid

i=average_sedimentation_rate_from_decomp_sed_thick_cm_per_ka_70_80.grd
j=`basename $i .grd`

grdsample $i -I1 -R-180/180/-70/80 -G${j}_1deg.grd
grd2xyz ${j}_1deg.grd >${j}.xyz

# raw sed thickness grid

# i=sedthick_world_v2.grd
i=Wobbe_Lindeque_Gohl_sedthick_world_v3_5min_GPC2014.grd
j=sedthick_1deg

if [ ! -s sedthick_1deg_masked.xyz ]; then
	echo "Converting $i to $j"
	grdsample $i -I1 -R-180/180/-70/80 -G${j}.grd
	grd2xyz ${j}.grd >${j}.xyz
	grd2xyz ${j}.grd | cut -f3 >${j}.z
fi

# convert sed rate grid

i=sed_rate_cm_per_ka.grd
j=`basename $i .grd`

if [ ! -s ${j}_1deg.z ]; then
	echo "Converting $i to $j.z"
	grdsample $i -I1 -R-180/180/-70/80 -G${j}_1deg.grd
	grd2xyz ${j}_1deg.grd | cut -f3 >${j}_1deg.z
fi


# convert river mouth mean distance grids to xyz and extract third colummn
# all output grids should be 1 deg pixel registered grids

for i in `ls *mean_distance*grd`

do
	j=`basename $i .grd`
	if [ ! -s ${j}_1deg.z ]; then
		echo "Converting $i to $j.z"
		grdsample -I1 -R-180/180/-70/80 -fg $i -G${j}_1deg.grd
		grd2xyz ${j}_1deg.grd | cut -f3 >${j}_1deg.z
	fi
done

