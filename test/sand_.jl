using ImplicitIntegration
using CartesianGeometry
using StaticArrays
const T=Float64
using Plots
Plots.default(show=true)
using Interpolations
using LinearAlgebra

# Parameters
nx, ny = 80, 80
x0, y0 = 0.0, 0.0
Lx, Ly = 1.0, 1.0
Δx, Δy = Lx/nx, Ly/ny

mesh = (collect(x0:Δx:Lx), collect(y0:Δy:Ly), [0.0, 0.01])
mesh_center = (0.5Δx:Δx:Lx-0.5Δx, 0.5Δy:Δy:Ly-0.5Δy)
println("Mesh Size: ", length(mesh[1]), " x ", length(mesh[2]))
println("Mesh Center Size: ", length(mesh_center[1]), " x ", length(mesh_center[2]))
Φ(x,y,t) = y - 0.2 - 0.05*sin(3π*x) + 0.1*sin(2π*t)

V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero)
A = CartesianGeometry.integrate(Tuple{1}, Φ, mesh, T, zero)
W = CartesianGeometry.integrate(Tuple{0}, Φ, mesh, T, zero, bary)
B = CartesianGeometry.integrate(Tuple{1}, Φ, mesh, T, zero, bary)

# Plot
Vₙ₊₁ = A[3][1:(nx+1)*(ny+1)]
Vₙ = A[3][(nx+1)*(ny+1)+1:end]
Vₙ = reshape(Vₙ, (nx+1, ny+1))'
Vₙ₊₁ = reshape(Vₙ₊₁, (nx+1, ny+1))'
heatmap(mesh[1], mesh[2], Vₙ, aspect_ratio=1, c=:viridis, color=:grays, grid=false, cbar=true)
readline()

# Calculate Height in each column
Hₙ = sum(Vₙ, dims=1)
Hₙ₊₁ = sum(Vₙ₊₁, dims=1)

println("Size of Height: ", size(Hₙ))
println("Height: ", Hₙ)

p1 = plot(mesh_center[1],Hₙ[:],seriestype=:scatter)
p2 = plot(mesh_center[1],Hₙ₊₁[:],seriestype=:scatter)
plot(p1,p2)

readline()

# Compute the interface position
sn = y0 .+ Hₙ ./ Δx
sn1 = y0 .+ Hₙ₊₁ ./ Δx

println("Size of s: ", size(sn))
println("s: ", sn)

p1 = plot(mesh_center[1],sn[:],seriestype=:scatter)
p2 = plot(mesh_center[1],sn1[:],seriestype=:scatter)
plot(p1,p2)

# Interface Term flux
q = 0.001*ones(size(Hₙ))

plot(mesh_center[1],q[:],seriestype=:scatter)

# Increase the Height 
new_Hn = Hₙ .+ q
new_Hn1 = Hₙ₊₁ .+ q

p1= plot(mesh_center[1],new_Hn[:],seriestype=:scatter)
p2= plot(mesh_center[1],new_Hn1[:],seriestype=:scatter)
plot(p1,p2)

# Compute the new interface position
sn = y0 .+ new_Hn ./ Δx
sn1 = y0 .+ new_Hn1 ./ Δx

p1= plot(mesh_center[1],sn[:],seriestype=:scatter)
p2= plot(mesh_center[1],sn1[:],seriestype=:scatter)
plot(p1,p2)

# Compute the space time level set
tₙ = 0.0
tₙ₊₁ = 0.01
Δt = tₙ₊₁ - tₙ
sₙ = linear_interpolation((mesh[1]), vec(sn), extrapolation_bc=Flat())
sₙ₊₁ = linear_interpolation((mesh[1]), vec(sn1), extrapolation_bc=Flat())

body = (xx,yy,tt)->(xx - (sₙ(yy)*(tₙ₊₁ - tt)/Δt + sₙ₊₁(yy)*(tt - tₙ)/Δt))

stmesh = (mesh[1], mesh[2], [tₙ, tₙ₊₁])

V, bary, interface_length, cell_type = CartesianGeometry.integrate(Tuple{0}, body, stmesh, T, zero)
