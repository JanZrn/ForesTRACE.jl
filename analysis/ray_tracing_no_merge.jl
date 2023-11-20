using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using CSV
using Statistics

### My modules
include("Voxels_Mod.jl")

using .VoxelsMod

############################################################
# Voxelisation of the environment and Analysis of the rays #
############################################################

##############################
# 1. Create a vector of rays #
##############################

### load merged ray data frame from "clip_and_merge_rays.jl"
df_merged_rays = DataFrame(CSV.File("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/dud/2023-10-13_09-46-29_dud_merged_pine_B.csv"))

### load single tree - for voxel space extent
header, points = LazIO.load("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/dud/2023-10-13_09-46-29_pine_B.laz")

### create a sample of 400_000 random rays
df_sampled_rays = df_merged_rays[sample(1:nrow(df_merged_rays), 400_000, replace=false), :]

### load/save sample df - if memory saving is needed
#CSV.write("source_file.csv", df_sampled_rays)
#df_sampled_rays = DataFrame(CSV.File("source_file.csv"))

### Creates vectors of sampled ray coordinates
rx0, ry0, rz0, rx1, ry1, rz1 = (df_sampled_rays.x), (df_sampled_rays.y), (df_sampled_rays.z), (df_sampled_rays.x_1), (df_sampled_rays.y_1), (df_sampled_rays.z_1)

### Create a vector of segments (rays) to be used in the analysis - work with a vector and not with a DataFrame!
ray_vector=[Segment((rx0[i], ry0[i], rz0[i]),(rx1[i], ry1[i], rz1[i])) for i in eachindex(rx0)]
stop_vector=[ray_vector[i].vertices[2] for i in eachindex(ray_vector)]


###########################################################
# 2. Create the voxel environment and perform ray tracing #
###########################################################

### Choose the size of the plot to be analysed (now 18x18xZ m with a 5 m buffer in +- X, Y)
### Outputs the boundaries of a given plot cloud in three dimensions
### End pts have a +- 5m buffer, but the voxel space is just the tree itself
extent_x, extent_y, extent_z = (header.x_min, header.x_max), (header.y_min, header.y_max), (header.z_min, header.z_max)

### Creates an indexable 3D grid of Voxels - Environment
### With defined dimensions and 0, 0 in the pass and stop value
### Choose desired Voxel size
env = create_voxels(extent_x,extent_y,extent_z, 0.5)

### Executes the ray tracin analysis and creates a DataFrame from it
df_raytrace = DataFrame(raytrace!(env, ray_vector, stop_vector))

### Here the analysed voxel space can be saved or the analysis can continue?
#CSV.write("source_file.csv", df_raytrace)

##################################################
# Computing metrics from the ray tracing results #
##################################################

### 1) Quantifying of openness
### Openness - emptiness of the space - if nothing passed or stopped in a voxel, it will be 0 and will not show in our voxel space (alpha = 0) etc
df_raytrace.openness_NaN = df_raytrace.pass ./ (df_raytrace.pass .+ df_raytrace.stop)
df_raytrace.openness .= [isnan(op) ? 0.0 : op for op in df_raytrace.openness_NaN]

### to quantify how open the scanned area is, I will simply divide the number of voxels with openness (?0.95?-1.0) by the number of analysed voxels 
### define the range, when we see the voxel as "open space"
lower_range = 0.95
upper_range = 1.00

### Extracts only the "open" voxels
df_raytrace_open = filter(row -> row.openness >= lower_range && row.openness <= upper_range, df_raytrace) #extracts voxels with openness > 0.95, could be useful later
### How much of the voxel space is open, in %
raytrace_openness = (nrow(df_raytrace_open) / size(df_raytrace,1)) * 100

### 2) Focus/user bias 
### Tells us if the space was scanned equally or if some parts were scanned more and other less
### Standard deviation and variability within .focus of the voxel space
df_raytrace.focus = ((df_raytrace.pass) ./ sum(df_raytrace.pass))

maximum(df_raytrace.focus)
std(df_raytrace.focus)
var(df_raytrace.focus)


### 3) Occlusion rate
### We can quantify how much of the scanned volume has been occluded
### If the space is occluded (value of .pass == 0), the value is 1, vice versa
df_raytrace.occlusion = float.(iszero.(df_raytrace.pass))

df_raytrace_occluded = filter(row -> row.occlusion == 1.0, df_raytrace)

occ_rate = (nrow(df_raytrace_occluded) / nrow(df_raytrace)) * 100


### Save the final DataFrame as .csv
CSV.write("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/dud/2023-10-13_09-46-29_raytraced_0.5v_pine_B.csv", df_raytrace)
###