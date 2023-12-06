using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using CSV
using Statistics
#using MeshViz
#using GLMakie

using BenchmarkTools


### My modules
include("PtCloudManipulation.jl")
include("VoxMod.jl")

using .PtCloudManpMod
using .VoxelsMod

### Load data
### DataFrame with raytraced voxels
df_raytrace = DataFrame(CSV.File("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_raytraced_0.25v_pine_A.csv"))

### Point cloud of the analysed tree
header, points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_pine_A.laz")

### Ground point cloud of the forest plot
ground_header, ground_points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_ground.laz")

### Recreate voxel space 
### Change this according to the voxel environment
voxels=create_voxels((header.x_min, header.x_max), (header.y_min, header.y_max), (header.z_min, header.z_max), 0.25) 

### Changes the string .poly into voxel .poly in the DataFrame
df_raytrace.poly = [voxels[i].poly for i in 1:length(voxels)]


### Function that goes through all the columns in the given environment and filters all voxels laying under the presumed ground
### It is fed a dataframe of the analysed environment and returns the filtered dataframe without voxels under ground
### Arguments are: DataFrame, voxel resolution and ground points with its header

### So far the fastest iteration, but I feel like it can still be optimized
### I played with broadcasting etc, but this was still the fastest

aboveground_vox = filter_underground_occ(df_raytrace, 0.25, ground_points, ground_header)

aboveground_vox.middles_x, aboveground_vox.middles_y, aboveground_vox.middles_z = get_middles(aboveground_vox)

gd = groupby(aboveground_vox, [:middles_x, :middles_y])

occlusion_rate(aboveground_vox)

maximum([top_vox(gd[g], 0.95) for g in 1:length(gd)])