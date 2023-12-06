module VoxelsMod

#################################################################
# Module that keeps all the voxel related functions and structs #
#################################################################

# Voxel is a 3D cube in space, it has a given side length 
# and position in space - which are contained in the .poly part.
# The two other struct info contain a number of passes - how many 
# rays passed through the voxel. And a number of stops - how many 
# rays ended in the particular voxel.

using Meshes
using TerminalLoggers
using ProgressLogging
using DataFrames
using FileIO, LasIO, LazIO
using Statistics
#using MeshViz
#using GLMakie

export Voxel, create_voxels, intersects, ray_voxel_intersect!, stop_voxel_intersect!, raytrace!
export get_middles, filter_underground_occ, occlusion_rate, top_vox
export voxel_viz_solids, voxel_viz_occlusion, voxel_viz_focus, voxel_viz_openness

### My modules
include("PtCloudManipulation.jl")
using .PtCloudManpMod


### The new struct Voxel
mutable struct Voxel 
    poly::Meshes.Box
    pass::Int32 # how many rays passed through the voxel
    stop::Int32 # how many rays ended in this voxel
end

### Creates a space of empty voxels - environment
function create_voxels(extent_x,extent_y,extent_z,voxel_side_length::Float64) 
    x0=ceil.(extent_x[1]) -1 
    y0=ceil.(extent_y[1]) -1
    z0=ceil.(extent_z[1]) -1 #-1 so it includes the edge laser-points
    xmax=ceil.(extent_x[2])
    ymax=ceil.(extent_y[2])
    zmax=ceil.(extent_z[2])
    grid=CartesianGrid((x0,y0,z0),(xmax,ymax,zmax),(voxel_side_length,voxel_side_length,voxel_side_length))

    [Voxel(Meshes.Box((grid[i].vertices[1]),(grid[i].vertices[7])), 0 , 0 ) for i in 1:length(grid)]
end

### A new iteration of the hasintersect function, that works with our Voxel structure
### If a Voxel is penetrated by the ray, does it also end there?

### Intersects voxels with rays
function intersects(box::Meshes.Box, seg::Meshes.Segment) # voxel and ray
    seg_start, seg_end = coordinates.(seg.vertices)
    box_min, box_max = coordinates(box.min), coordinates(box.max)

    seg_length = seg_end .- seg_start
    seg_start_to_box_min = box_min .- seg_start
    seg_start_to_box_max = box_max .- seg_start

    tnear = floatmin(Float64)
    tfar = floatmax(Float64)

    for axis in eachindex(seg_length)
        if seg_length[axis] == 0 # parallel
            if (seg_start_to_box_min[axis] > 0) || (seg_start_to_box_max[axis] < 0)
                    return false # segment is not between planes
            end
        else
            t1 = seg_start_to_box_min[axis] / seg_length[axis]
            t2 = seg_start_to_box_max[axis] / seg_length[axis]

            tmin = minimum([t1, t2])
            tmax = maximum([t1, t2])

            if tmin > tnear   tnear = tmin end
            if tmax < tfar    tfar = tmax end
            if (tnear > tfar) || (tfar < 0.0)
                return false
            end
        end
    end
    if (tnear >= 0) && (tnear <= 1)  return true end
    if (tfar >= 0) && (tfar <= 1)  return true end
    return false
end

### Intersects voxels with end points
function intersects(box::Meshes.Box{3, Float64}, point::Meshes.Point) # voxel and the end point
    x, y, z = coordinates(point)
    xmax, ymax, zmax = coordinates(box.max)
    xmin, ymin, zmin = coordinates(box.min)
    ((x < xmax) && (x > xmin)) && ((y < ymax) && (y > ymin)) && ((z < zmax) && (z > zmin))
end

### Cycles the intersect function with ::Segment through the voxel space
function ray_voxel_intersect!(voxel::Voxel,ray::Meshes.Segment) # adds the .pass value into the Voxel struct
    didPass = intersects(voxel.poly, ray) # did a given ray intersect a given voxel? (0/1)
    voxel.pass += didPass # adds the (0/1) to .pass in our Voxel structure
end

### Cycles the intersect function with ::Point through the voxel space
function stop_voxel_intersect!(voxel::Voxel, stop::Meshes.Point) # adds the .stop value into the Voxel struct
    didStop = intersects(voxel.poly, stop) # if it ends there = 1
    voxel.stop += didStop # adds (0/1) to .stop in our Voxel structure
end

### Function that loops through all the given voxels and rays+end-points 
### Single-thread with a progress bar
function raytrace!(voxel_space::Vector{Main.VoxelsMod.Voxel}, ray_vector::Vector{Segment{3, Float64, Vector{Meshes.Point3}}}, stop_vector::Vector{Meshes.Point3})
    @progress for i in eachindex(voxel_space)
        ray_voxel_intersect!.((voxel_space[i], ), ray_vector) 
        stop_voxel_intersect!.((voxel_space[i], ), stop_vector) 
    end
    return voxel_space # returns the modified voxel space
end

### Adds middles of voxels into the DataFrame
### Syntax = df.middles_x, df.middles_y, df.middles_z = get_middles(df)
function get_middles(df::DataFrame)
    df.middles_x .= [(coordinates(df.poly[i](0))[1] + ((coordinates(df.poly[i](0))[1] - coordinates(df.poly[i](1))[1])/2)) for i in 1:size(df, 1)]
    df.middles_y .= [(coordinates(df.poly[i](0))[2] + ((coordinates(df.poly[i](0))[2] - coordinates(df.poly[i](1))[2])/2)) for i in 1:size(df, 1)]
    df.middles_z .= [(coordinates(df.poly[i](0))[3] + ((coordinates(df.poly[i](0))[3] - coordinates(df.poly[i](1))[3])/2)) for i in 1:size(df, 1)]
    return df.middles_x, df.middles_y, df.middles_z
end

### Function that goes through all the columns in the given environment and filters all voxels laying under the presumed ground
### It is fed a dataframe of the analysed environment and returns the filtered dataframe without voxels under ground

### So far the fastest iteration, but I feel like it can still be optimized
### I played with broadcasting etc, but this was still the fastest

function filter_underground_occ(df::DataFrame, voxel_resolution, ground_points, ground_header) # Arguments are: DataFrame, voxel resolution and ground points with its header

    ### Create individual columns of set X and Y and changing Z from the DataFrame
    ### Extracting the coordinates of voxel middles
    df.middles_x, df.middles_y, df.middles_z = get_middles(df)
    
    ### Creates columns (groups) with same X and Y coordinates
        gdf = groupby(df, [:middles_x, :middles_y] )
        x = [gdf[i].middles_x[1] for i in eachindex(gdf)] # extracts x coords of each column
        y = [gdf[i].middles_y[1] for i in eachindex(gdf)] # extracts y coords of each column
        kernel = [Meshes.Box((x[i] - (1.5*voxel_resolution), y[i] - (1.5*voxel_resolution)), (x[i] + (1.5*voxel_resolution), y[i] + (1.5*voxel_resolution))) for i in eachindex(x)] # creates a 3x3 (*vox size) kernel for each column
        
        ground_pts = [filter_pixel!(ground_points, ground_header, kernel[i]) for i in eachindex(kernel)] # filters ground points for each column
        ground_z = [mean(zcoord.(ground_pts[i], (ground_header,))) for i in eachindex(ground_pts)] # gets a mean ground value for each column
    
        grnd_coord = DataFrame(:x => x, :y => y, :z => ground_z) # creates a Dataframe for each column
        filtered_groups = [filter(row -> row.middles_z >= ground_z[i], gdf[(middles_x = grnd_coord.x[i], middles_y = grnd_coord.y[i])]) for i in eachindex(ground_z)] # deletes voxels below ground value
        
        voxels_no_ground_occlusion = vcat(filtered_groups...) # merges all into one DataFrame
        
return voxels_no_ground_occlusion
end

### Occlusion rate function
function occlusion_rate(df::DataFrame)
    occ_rate::Float64 = (nrow(filter(row -> row.occlusion == 1.0, df)) / nrow(df)) * 100
return occ_rate
end

### For G-T method of acquiring the highest non-open voxel
### Must be a DataFrame
### Cycles through all the colummns and finds the highest non-open voxel (threshold < 0,95)
function top_vox(gdf, threshold::Float64)
    topvox = 0.0
    x = 0.0
    y = 0.0

for i in 1:size(gdf, 1)
    if (0.0 < gdf.openness[i] < threshold) & (gdf.stop[i] > 1) & (gdf.occlusion[i] == 0.0) # if the voxel is considered open, has at least two stop points and non occluded
        topvox = gdf.middles_z[i]
        x = gdf.middles_x[i]
        y = gdf.middles_y[i]
    else 
        0
    end
end
return topvox, x, y #returns (Z, X, Y)
end


### Visualization
### Visualization of openness
### If solid = true, will show only the solids with their openness, if false, will show evetrything
function voxel_viz_openness(voxel::DataFrame, solid::Bool, threshold)
    p = Scene()

    if solid == true
        voxel! = filter(row -> row.openness < threshold && row.occlusion == 0.0, voxel)
    p = viz(voxel!.poly, color = voxel!.openness, alpha = (voxel!.openness))
    else
    p = viz(voxel.poly, color = voxel.openness, alpha = (voxel.openness))
    end
    return p
end

### Visualization of focus/user bias
function voxel_viz_focus(voxel::DataFrame)
    p = Scene()
    p = viz(voxel.poly, color = voxel.focus , alpha = (voxel.focus))
    return p
end

### Visualization of occlusion
function voxel_viz_occlusion(voxel::DataFrame, clr)
    voxel! = filter(row -> row.occlusion == 1.0, voxel)
    p = Scene()
    p = viz(voxel!.poly, color = clr, alpha = 0.5)
    return p
end

### Visualization of voxels with openness lower than a threshold for openness
function voxel_viz_solids(voxel::DataFrame, clr, threshold)
    voxel! = filter(row -> row.openness < threshold && row.occlusion == 0.0, voxel)
    p = Scene()
    p = viz(voxel!.poly, color = clr, alpha = 0.5)
    return p
end

end # end of module