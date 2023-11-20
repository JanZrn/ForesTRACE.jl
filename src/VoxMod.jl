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

export Voxel, create_voxels, ray_voxel_intersect!, stop_voxel_intersect!, raytrace!, top_vox

### The new struct Voxel
mutable struct Voxel 
    poly::Meshes.Box
    pass::Int32 # how many rays passed through the voxel
    stop::Int32 # how many rays ended in this voxel
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
    if (0.0 < gdf.openness[i] < threshold) & (gdf.occlusion[i] == 0.0) # if the voxel is considered open and non occluded
        topvox = gdf.middles_z[i]
        x = gdf.middles_x[i]
        y = gdf.middles_y[i]
    else 
        0
    end
end
return topvox, x, y #returns (Z, X, Y)
end

end # end of module