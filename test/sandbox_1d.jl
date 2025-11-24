using CartesianGeometry
using Test
const T=Float64

universe = (-1:11,)
node = (1:9,)

# define mesh
xyz = collocated.(identity, universe, node)

@assert all(@. isequal(length(xyz), length(universe)))

# define level set
const R = 0.2
const a = 0.5
nx = 10
dx = 1.0/nx
x0 = 0.0

mesh = ([x0 + i*dx for i in 0:nx-1],)

levelset = HyperSphere(R, (a,))
levelset = (x, _=0) -> (x - a)

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, mesh, T, zero)
As = integrate(Tuple{1}, levelset, mesh, T, zero)

Ws = integrate(Tuple{0}, levelset, mesh, T, zero, bary)
Bs = integrate(Tuple{1}, levelset, mesh, T, zero, bary)

@show interface_length
@show length(interface_length)
# Test ImplicitIntegration
nx = 10
dx = 1.0/nx
x0 = 0.0
mesh = ([x0 + i*dx for i in 0:nx-1],)

levelset = (x) -> (x[1] - 0.5)

V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh, levelset)

@show Γ1
@show length(Γ1)

@test cell_types1 == cell_types