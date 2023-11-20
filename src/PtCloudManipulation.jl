module PtCloudManpMod

######################################################################
# Module created to keep all the new functions dealing with 
# pt clouds and manipulating with their data. 

using DelimitedFiles
using FileIO, LasIO, LazIO
using Meshes

export within, filter_pixel!

### Tells what laz points are inside of a given geometry
function within(point::LasPoint, header::LasHeader, pixel::Meshes.Box{2, Float64})
    x, y = xcoord(point, header), ycoord(point, header)
    xmax, ymax = coordinates(pixel.max)
    xmin, ymin = coordinates(pixel.min)
    ((x < xmax) && (x > xmin)) && ((y < ymax) && (y > ymin))
end

### Returns .laz points that are in a given geometry
function filter_pixel!(points::Vector, header::LasHeader, pixel::Meshes.Box{2, Float64})
    points[within.(points, (header, ), (pixel, ))]
end




end # end of module