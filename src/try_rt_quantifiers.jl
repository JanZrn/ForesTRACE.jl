using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using CSV
using Statistics

### My modules
include("PtCloudManipulation.jl")
include("VoxMod.jl")

using .PtCloudManpMod
using .VoxelsMod

### Load data
### DataFrame with raytraced voxels
df_raytrace = DataFrame(CSV.File("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_raytraced_0.25v_pine_B.csv"))

### Point cloud of the analysed tree
header, points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_pine_B.laz")

### Ground point cloud of the forest plot
ground_header, ground_points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_ground.laz")

voxels=create_voxels((header.x_min, header.x_max), (header.y_min, header.y_max), (header.z_min, header.z_max), 0.25) 

df_raytrace.poly = [voxels[i].poly for i in 1:length(voxels)]

df_raytrace_n = select(df_raytrace,Not([:occlusion,:focus, :openness, :openness2]))

rt_quantifiers(df_raytrace_n)

filter_underground_occ(df_raytrace_n, 0.25, ground_points, ground_header)

get_topvox(df_raytrace_n, 0.95)

gt_height(df_raytrace_n, 0.95, ground_points, ground_header)

println(df_raytrace_n[69000:69005, :])