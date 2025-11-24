using CartesianGeometry
using Test
const T=Float64
using StaticArrays

x = collect(0.0:0.05:1.0)
y = collect(0.0:0.05:1.0)
z = [0.0, 1.0]
xyz = (x, y, z)

# define level set
# Define initial and final circle parameters
const R_init = 0.1    # Initial radius
const R_final = 0.10000000001   # Final radius 
const a, b = 0.5, 0.5 # Circle center

levelset1 = (x, y, z) -> (x-a)^2 + (y-b)^2 - R_init^2
levelset2 = (x, y, z) -> (x-a)^2 + (y-b)^2 - R_final^2

tn,tn1 = 0.0, 1.0
dt = tn1 - tn
levelset = (x, y, t) -> (t-tn)/dt * levelset2(x, y, 0) + (1 - (t-tn)/dt) * levelset1(x, y, 0)

V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
As = integrate(Tuple{1}, levelset, xyz, T, nan)

Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)
