
#########################################################################
# Calculating the height of a given tree after the ray tracing analysis #
#########################################################################

using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using CSV
using Statistics

### My modules
include("PtCloudManipulation_Mod.jl")
include("Voxels_Mod.jl")

using .PtCloudManpMod
using .VoxelsMod

### Load data
### DataFrame with raytraced voxels
df_raytrace = DataFrame(CSV.File("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_raytraced_0.5v_pine_B.csv"))

### Point cloud of the analysed tree
header, points = LazIO.load("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_pine_B.laz")

### Ground point cloud of the forest plot
ground_header, ground_points = LazIO.load("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_ground.laz")

### Recreate voxel space 
### Change this according to the voxel environment
voxels=create_voxels((header.x_min, header.x_max), (header.y_min, header.y_max), (header.z_min, header.z_max), 0.5) 

### Changes the string .poly into voxel .poly in the DataFrame
df_raytrace.poly = [voxels[i].poly for i in 1:length(voxels)]


### Create individual columns of set X and Y and changing Z from the DataFrame
### Extracting the coordinates of voxel middles
df_raytrace.middles_x .= [(coordinates(df_raytrace.poly[i](0))[1] + ((coordinates(df_raytrace.poly[i](0))[1] - coordinates(df_raytrace.poly[i](1))[1])/2)) for i in 1:size(df_raytrace, 1)]
df_raytrace.middles_y .= [(coordinates(df_raytrace.poly[i](0))[2] + ((coordinates(df_raytrace.poly[i](0))[2] - coordinates(df_raytrace.poly[i](1))[2])/2)) for i in 1:size(df_raytrace, 1)]
df_raytrace.middles_z .= [(coordinates(df_raytrace.poly[i](0))[3] + ((coordinates(df_raytrace.poly[i](0))[3] - coordinates(df_raytrace.poly[i](1))[3])/2)) for i in 1:size(df_raytrace, 1)]


### Creates columns (groups) with same X and Y coordinates
gd = groupby(df_raytrace, [:middles_x, :middles_y] )

####################################################
# Ground - Topvox method for measuring tree height #
####################################################

### Topvox - the highest non-open voxel (openness < 0.95)
### Tells us what middle-coordinates the highest non-open voxel has (z, x, y)
topvox_info = maximum([top_vox(gd[g], 0.95) for g in 1:length(gd)])

### Finds Ground value (Z) bellow the Topvox
### Create a 1x1 m buffer on the ground below the Topvox
ground_buffer = Meshes.Box((topvox_info[2] - 1, topvox_info[3] - 1), (topvox_info[2] + 1, topvox_info[3] + 1))

### Clip the ground points
ground_points, ground_header = filter_pixel!(ground_points, ground_header, ground_buffer)

### Get the mean Z value of ground
ground_z = mean(zcoord(ground_points, header))


### The tree height is equal to the distance between ground and the center of our topvox (+- voxel side length / 2)
G_T_height = (topvox_info[1] + abs(ground_z))

