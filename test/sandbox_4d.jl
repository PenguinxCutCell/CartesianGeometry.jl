using CartesianGeometry
using Test
const T = Float64

universe = (-1:11, -1:11, -1:11, -1:11)
node = (1:9, 1:9, 1:9, 1:9)

# define mesh
xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.25
const a, b, c, d = 0.5, 0.5, 0.5, 0.5

levelset = HyperSphere(R, (a, b, c, d))

V, bary, interface_volume, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)

@show cell_types
@show interface_volume