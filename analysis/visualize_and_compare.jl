using FileIO, LasIO, LazIO
using Meshes
using DataFrames
using MeshViz
using GLMakie
using CSV
using Colors, ColorSchemes

### My modules
include("Voxels_Mod.jl")

using .VoxelsMod

### Since I couldn't find a way how to properly save the Voxel DataFrame, so the .poly is type Meshes.Box and not a string
### I decided to recreate the voxel space and pass it into the DataFrames

### load the DataFrames to be compared
df_raytrace = DataFrame(CSV.File("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_raytraced_0.1v_pine_B.csv"))

### Point cloud of the analysed tree
header, points = LazIO.load("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/tops/2023-10-13_09-44-40_pine_B.laz")

### Recreate voxel space 
### Change this according to the voxel environment
voxels=create_voxels((header.x_min, header.x_max), (header.y_min, header.y_max), (header.z_min, header.z_max), 0.1) 

### Changes the string .poly into voxel .poly in the DataFrame
df_raytrace.poly = [voxels[i].poly for i in 1:length(voxels)]


### Visualization of openness
function voxel_viz_openness(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = voxel.openness, alpha = (voxel.openness))
    return p
end

### Visualization of focus/user bias
function voxel_viz_focus(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = voxel.focus , alpha = (voxel.focus))
    return p
end

### Visualization of occlusion
function voxel_viz_occlusion(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = :red, alpha = 0.5)
    return p
end

### Visualization of .stop
function voxel_viz_stop(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = :green, alpha = 0.5)
    return p
end

function voxel_viz_grey(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = :grey, alpha = (voxel.openness))
    return p
end

### visualization

### 1) Quantifying of openness
### to quantify how open the scanned area is, I will simply divide the number of voxels with openness (?0.95?-1.0) by the number of analysed voxels 

### define the range, when we see the voxel as "open space"
lower_range = 0.95
upper_range = 1.00

### Extracts only the "open" voxels
df_raytrace_open = filter(row -> row.openness2 >= lower_range && row.openness2 <= upper_range, df_raytrace) #extracts voxels with openness > 0.95, could be useful later
### How much of the voxel space is open, in %
raytrace_openness = (nrow(df_raytrace_open) / length(environment)) * 100

openness = voxel_viz_openness(df_raytrace)
save("D:/jazr0001/Ray_Tracing_MasterThesis/data/plots/swe/2022-10-06_09-57-12/1x1x1_400.000r_2022-10-06_09-57-12_openness.png", openness)

### 2) Focus/user bias - a square root of the number is taken to make the difference more visible
focus = voxel_viz_focus(df_raytrace)
save("D:/jazr0001/Ray_Tracing_MasterThesis/data/plots/swe/2022-10-06_09-57-12/1x1x1_400.000r_2022-10-06_09-57-12_focus.png", focus)

### 3) Occlusion rate
### If the space is occluded (value of .pass == 0), the value is 1, vice versa
df_raytrace_occluded = filter(row -> row.occlusion == 1.0, df_raytrace)
occ_rate = (nrow(df_raytrace_occluded) / nrow(df_raytrace)) * 100

occlusion = voxel_viz_occlusion(df_raytrace_occluded)
save("D:/jazr0001/Ray_Tracing_MasterThesis/data/plots/swe/2022-10-06_09-57-12/1x1x1_400.000r_2022-10-06_09-57-12_occlusion.png", occlusion)


### 4) Plot Visualization
### Just to show the plot - showing the .stop
df_raytrace_stop = filter(row -> row.stop > 0, df_raytrace)

plot_stop = voxel_viz_stop(df_raytrace_stop)
save("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/pine_B_stop_control.png", plot_stop)

### Non-open
df_raytrace_nonopen = filter(row -> row.openness < 0.95 && row.occlusion == 0.0, df_raytrace)
plot_nonopen = voxel_viz_openness(df_raytrace_nonopen)
save("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/pine_B_open_control.png", plot_nonopen)

### voxel_viz_grey

grey = voxel_viz_grey(df_raytrace_nonopen)
save("D:/jazr0001/Ray_Tracing_MasterThesis/95th_percentile/G-T_grey.png", grey)