using CartesianGeometry
using Test
const T=Float64

universe = (-1:11, -1:19)
#
#halos = ((-1:10, -1:18),
#         (0:10, 0:18))

node = (1:9, 1:17)
#center = (1:8, 1:16)

# define mesh
xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.25
const a, b = 0.5, 0.5

levelset = HyperSphere(R, (a, b))

#=
vol, bary = integrate(Tuple{0}, levelset, xyz, T, zero)

@assert isequal(length(vol), prod(length.(universe)))
@assert isequal(length(bary), prod(length.(universe)))

@assert isapprox(sum(vol), π * R ^ 2)
@assert isapprox(sum(first.(bary) .* vol), π * R ^ 2 * a)
@assert isapprox(sum(last.(bary) .* vol), π * R ^ 2 * b)

surf = integrate(Tuple{1}, levelset, xyz, T, zero)

@assert all(isequal.(length.(surf), prod(length.(universe))))
=#

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@test size(V) == size(bary) == size(interface_length)

@show V
@show bary
@show interface_length

# Calculate Perimeter of the circle
# Remove NaN values
interface_length = interface_length[.!isnan.(interface_length)]
println(sum(interface_length))
println(2 * π * R)

@assert isapprox(sum(interface_length), 2 * π * R, atol=1e-7)

# Cell types
@show cell_types

# Verify that indices with cell_types == -1 have corresponding interface_length values
indices = findall(cell_types .== -1)
indices_inter = findall(interface_length .> 0)
@show length(indices)
@show length(indices_inter)

# Print me only th interface_length[indices] values
@show interface_length[indices]
@show interface_length[indices_inter]

# second-kind moments
#tmp = integrate(Tuple{0}, levelset, xyz, bary)

#=
using StaticArrays
using CartesianArrays
using CartesianGeometry

import Base: OneTo

const T = Float64

#=== Domains ===#
universe = (-1:11, -1:19)

halos = ((-1:10, -1:18),
         (0:10, 0:18))
node = (1:9, 1:17)
center = (1:8, 1:16)

#=== Mesh ===#
xyz = collocated.(Ref(identity), universe, node)

#=== Geometry ===#
levelset = HyperSphere(0.25, 0.5 .* one.(eltype.(xyz)))

#=== Capacities (1st kind) ===#
#--- Volume ---#
mom = CVector{SVector{length(halos[1])+1,T}}(undef, halos[1])
integrate!(mom, Tuple{0}, levelset, xyz, halos[1])

v = first.(mom)
bary = deleteat.(mom, 1)
mom = nothing

#--- Surface ---#
domains = ntuple(_ -> halos[1], length(halos[1]))
a = CVector{T}(undef, (halos[1]..., OneTo(length(halos[1]))))
integrate!(a, Tuple{1}, levelset, xyz, domains)

#=== Capacities (2nd kind) ===#
#--- Volume ---#
domains = ntuple(_ -> halos[2], length(halos[2]))
w = CVector{T}(undef, (halos[2]..., OneTo(length(halos[2]))))
integrate!(w, Tuple{0}, levelset, xyz, bary, domains)

#--- Surface ---#
b = CVector{T}(undef, (halos[1]..., OneTo(length(halos[1]))))
integrate!(b, Tuple{1}, levelset, xyz, bary, halos[1])

=#
