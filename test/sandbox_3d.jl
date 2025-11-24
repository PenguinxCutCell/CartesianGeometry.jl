using CartesianGeometry
using Test
const T=Float64

# Test 3D
universe = (-1:11, -1:11, -1:11)
node = (1:9, 1:9, 1:9)

xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.25
const a, b, c = 0.5, 0.5, 0.5

levelset = HyperSphere(R, (a, b, c))

V, bary, interface_area, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@test typeof(Bs) == Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

@test size(V) == size(bary) == size(interface_area)

@show interface_area

# Calculate Surface Area of the sphere
# Remove NaN values
interface_area = interface_area[.!isnan.(interface_area)]
println(sum(interface_area))
println(4 * π * R ^ 2)

@test isapprox(sum(interface_area), 4 * π * R ^ 2, atol=1e-2)

# Cell types
@show cell_types

# Verify that indices with cell_types == -1 have corresponding interface_area values
indices = findall(cell_types .== -1)
indices_inter = findall(interface_area .> 0)
@show length(indices)
@show length(indices_inter)

# Print me only th interface_area[indices] values
@show interface_area[indices]
@show interface_area[indices_inter]