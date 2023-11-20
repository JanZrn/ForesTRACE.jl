using DelimitedFiles
using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using CSV
using Dates
using SortMerge

### My modules
include("PtCloudManipulation_Mod.jl")
include("Voxels_Mod.jl")

using .PtCloudManpMod
using .VoxelsMod

#################################################
# Clip and Merge rays from the analysed ptcloud #
#################################################

### LOAD POINT CLOUD - ends
header, points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/bottoms/2023-10-13_09-43-09.laz") #load end position 

### Clipping of the point cloud into our desired size
### Size of the clipped ptcloud - Bounding Box
BB = Meshes.Box((-5.0, -5.0), (23.0, 23.0)) ### 5m buffer (the analysed plot is to be (0.0, 0.0), (18.0, 18.0))

## clip points into fpouints
fpoints = filter_pixel!(points, header, BB)

### LOAD TRAJECTORIES
data, header2 = readdlm("D:/jazr0001/Ray_Tracing_MasterThesis/data/swe/SE_160/2022-10-06_09-40-25_results/2022-10-06_09-40-25.gs-traj", ' ',header=true)


### Creating a DataFrame of rays in the ptcloud

# Convert the DateTime of end points
time_raw_ends=Dates.unix2datetime.(gps_time.(fpoints))

# get coordinates of ends
coordinates_ends=[[xcoord(fpoints[i], header ), ycoord(fpoints[i], header), zcoord(fpoints[i], header)] for i in 1:length(fpoints)]

# DataFrame of the end points df_ends
df_ends=DataFrame(world_time=time_raw_ends, coordinates=coordinates_ends)
df_ends.x = getindex.(df_ends.coordinates, 1)
df_ends.y = getindex.(df_ends.coordinates, 2)
df_ends.z = getindex.(df_ends.coordinates, 3)
select!(df_ends, Not(:coordinates))
df_ends 

# create DataFrame of trajectories
traj_raw=DataFrame(data, vec(header2))
# rename
traj_raw_rightname=rename(traj_raw, 1 => :pworld_time)
# extract and convert time into value
time_raw_traj=traj_raw_rightname.pworld_time
time_traj_date=unix2datetime.(time_raw_traj)

# add scaled time
df_traj=insertcols!(traj_raw_rightname, 2, :world_time => time_traj_date)
select!(df_traj, Not(:pworld_time, :q0, :q1, :q2, :q3, :r, :g, :b, :nx, :ny, :nz, :roll, :pitch, :yaw))



### save - if needed
CSV.write("df_traj_ray_tracing_2022-10-04_08-13-20.csv", df_traj)
CSV.write("df_ends_ray_tracing_2022-10-04_08-13-20.csv", df_ends)

### load - if needed
df_traj_load = DataFrame(CSV.File("df_traj_ray_tracing_2022-10-04_08-13-20.csv"))
df_ends_load = DataFrame(CSV.File("df_ends_ray_tracing_2022-10-04_08-13-20.csv"))



### Fuzzy merge, with a threshold (Millisecond(XXX))
### (x,y,z) = traj, (x_1,y_1,z_1) = ends 
k = sortmerge(df_traj.world_time, df_ends.world_time, Millisecond(50),
sd=(v1, v2, i1, i2, threshold) -> begin
diff = v1[i1] - v2[i2]
(abs(diff) > threshold)  &&  (return sign(diff))
return 0
end)

### merged DataFrame
df_merged_n = hcat(df_traj[k[1], :], df_ends[k[2], :], makeunique=true)

### Save the filtered(clipped) and merged DataFrame of rays
CSV.write("D:/jazr0001/Ray_Tracing_MasterThesis/data/swe/SE_160/2022-10-06_09-40-25_results/2022-10-06_09-40-25_df_merged_rayblock.csv", df_merged_n)      


############################################################
# Voxelisation of the environment and Analysis of the rays #
############################################################

##############################
# 1. Create a vector of rays #
##############################

### load merged ray data frame from "clip_and_merge_rays.jl"
df_merged_rays = DataFrame(CSV.File("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/bottoms/2023-10-13_09-43-09_bottoms_merged_pine_A.csv"))

### load single tree - for voxel space extent
header, points = LazIO.load("E:/Ray_Tracing/Ray_Tracing_MasterThesis/95th_percentile/bottoms/2023-10-13_09-43-09_pine_A.laz")

### create a sample of 400_000 random rays
df_sampled_rays = df_merged_rays[sample(1:nrow(df_merged_rays), 400_000, replace=false), :]

### load/save sample df - if memory saving is needed
#CSV.write("source_file.csv", df_sampled_rays)
#df_sampled_rays = DataFrame(CSV.File("source_file.csv"))

### Creates vectors of sampled ray coordinates
### ends
rx0 = (df_sampled_rays.x)
ry0 = (df_sampled_rays.y)
rz0 = (df_sampled_rays.z)
### traj
rx1 = (df_sampled_rays.x_1)
ry1 = (df_sampled_rays.y_1)
rz1 = (df_sampled_rays.z_1)

### Create a vector of segments (rays) to be used in the analysis - work with a vector and not with a DataFrame!
ray_vector=[Segment((rx0[i], ry0[i], rz0[i]),(rx1[i], ry1[i], rz1[i])) for i in eachindex(rx0)]
stop_vector=[ray_vector[i].vertices[2] for i in eachindex(ray_vector)]

###########################################################
# 2. Create the voxel environment and perform ray tracing #
###########################################################

### Choose the size of the plot to be analysed (now 18x18xZ m with a 5 m buffer in +- X, Y)
### Outputs the boundaries of a given plot cloud in three dimensions
extent_x=(header.x_min, header.x_max) #ptcloud.x has a +-5m buffer
extent_y=(header.y_min, header.y_max) #ptcloud.y has a +-5m buffer
extent_z=(header.z_min, header.z_max)


### Creates an indexable 3D grid of Voxels - Environment
### With defined dimensions and 0, 0 in the pass and stop value
### Choose desired Voxel size
env = create_voxels(extent_x,extent_y,extent_z, 5.0)

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
#CSV.write("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/static/2023-10-13_09-38-06_raytraced_0.25v_pine_A.csv", df_raytrace)
###
