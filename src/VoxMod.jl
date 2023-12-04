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
using MeshViz
using GLMakie

export Voxel, create_voxels, intersects, ray_voxel_intersect!, stop_voxel_intersect!, raytrace!, top_vox, voxel_viz_solids, voxel_viz_occlusion, voxel_viz_focus, voxel_viz_openness

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

### For G-T method of acquiring the highest non-open voxel
### Must be a GroupedDataFrame
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