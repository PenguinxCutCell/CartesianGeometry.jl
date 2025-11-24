using Test
using CartesianGeometry
using StaticArrays  # if required by CartesianGeometry

const T = Float64

@testset "1D Geometry Integration" begin
    # Setup mesh and level set
    R = 0.2
    a = 0.42
    nx = 10
    dx = 1.0/nx
    x0 = 0.0
    mesh = ([x0 + i*dx for i in 0:nx-1],)

    # Use a HyperSphere then override with a custom levelset function
    levelset = HyperSphere(R, (a,))
    levelset = (x, _=0) -> (x - a)

    # Compute moments using integrate for Tuple{0} and Tuple{1}
    V, bary, interface_length, cell_types = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero)
    As = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero)
    Ws = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero, bary)
    Bs = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero, bary)

    # Now test the implicit integration (new feature)
    mesh2 = ([x0 + i*dx for i in 0:nx-1],) # same mesh
    levelset2 = (x) -> (x[1] - a)  # simple 1D levelset function

    V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh2, levelset2)
    
    # Check that the results are the same
    @test length(V) == length(V1)
    @test length(cell_types) == length(cell_types1)
    @test length(interface_length) == length(Γ1)
    @test length(As) == length(A1)
    @test length(Ws) == length(W1)
    @test length(Bs) == length(B1)

    @test all(@. isapprox(V, V1, atol=1e-12))
    @test all(@. isapprox(cell_types, cell_types1, atol=1e-12))
    @test all(@. isapprox(interface_length, Γ1, atol=1e-12))
    @test all(@. isapprox(As, A1, atol=1e-12))  
    @test all(@. isapprox(Ws, W1, atol=1e-12))
    @test all(@. isapprox(Bs, B1, atol=1e-12))

end