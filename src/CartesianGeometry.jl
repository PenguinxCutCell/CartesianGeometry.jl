module CartesianGeometry

using Base.Cartesian
using StaticArrays
using TiledIteration
using Vofinit
using VofiJul
#using CartesianCore
#using CartesianArrays

using ImplicitIntegration

export nan
export HyperSphere
export collocated, staggered
export integrate, integrate!, get_cell_type
export implicit_integration

include("utils.jl")
include("zoo.jl")
include("mesh.jl")
include("vofinit.jl")
include("first.jl")
include("second.jl")
include("implicitint.jl")


const VOFMethod = Vofinit
const VOFIJulMethod = VofiJul
end