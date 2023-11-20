using DelimitedFiles
using FileIO, LasIO, LazIO
using Meshes
using Dates
using DataFrames
using CSV
using SortMerge

### My modules
include("PtCloudManipulation_Mod.jl")

using .PtCloudManpMod


### LOAD POINT CLOUD
header, points = LazIO.load("D:/jazr0001/Ray_Tracing_MasterThesis/data/swe/2022-10-04_08-23-59.laz") #load end position 

### Clipping of the point cloud into our desired size
function within(point::LasPoint, header::LasHeader, pixel::Box{2, Float64})
    x, y = xcoord(point, header), ycoord(point, header)
    xmax, ymax = coordinates(pixel.max)
    xmin, ymin = coordinates(pixel.min)
    ((x < xmax) && (x > xmin)) && ((y < ymax) && (y > ymin))
end

function filter_pixel(points::Vector, header::LasHeader, pixel::Box{2, Float64})
    points[within.(points, (header, ), (pixel, ))]
end

### size of the clipped ptcloud
BB = Meshes.Box((-5.0, -5.0), (23.0, 23.0)) ### 5m buffer (the analysed plot is to be (0.0, 0.0), (18.0, 18.0))

## clip points into fpouints
fpoints = filter_pixel!(points, header, BB)

### LOAD TRAJECTORIES
data, header2 = readdlm("D:/jazr0001/Ray_Tracing_MasterThesis/data/swe/2022-10-04_08-23-59.gs-traj", ' ',header=true)


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



### fuzzy merge, with a threshold (Millisecond(XXX))
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
CSV.write("D:/jazr0001/Ray_Tracing_MasterThesis/data/swe/2022-10-04_08-23-59_df_merged_rayblock.csv", df_merged_n)      


