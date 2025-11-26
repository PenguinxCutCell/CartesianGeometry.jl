module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit
using VofiJul
#using CartesianCore
#using CartesianArrays

export nan
export HyperSphere
export collocated, staggered
export integrate, integrate!, integrate_centroid, get_cell_type

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("vofinit.jl")
include("first.jl")
include("second.jl")


const VOFMethod = Vofinit
const VOFIJulMethod = VofiJul
end
