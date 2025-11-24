using Test
using CartesianGeometry
using StaticArrays

const T = Float64

@testset "3D Geometry Integration" begin
    # Setup mesh and level set
    R = 0.21
    a = 0.5
    nx, ny, nz = 8, 8, 8
    dx, dy, dz = 1.0/nx, 1.0/ny, 1.0/nz
    x0, y0, z0 = 0.0, 0.0, 0.0
    mesh = (
        [x0 + i*dx for i in 0:nx-1], 
        [y0 + i*dy for i in 0:ny-1],
        [z0 + i*dz for i in 0:nz-1]
    )

    # Use a 3D sphere as the level set function
    levelset = HyperSphere(R, (a, a, a))
    levelset = (x, y, z, _=0) -> (x - a)^2 + (y - a)^2 + (z - a)^2 - R^2

    # Compute moments using integrate for Tuple{0} and Tuple{1}
    V, bary, interface_area, cell_types = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero)
    As = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero)
    Ws = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero, bary)
    Bs = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero, bary)

    # Now test the implicit integration (new feature)
    mesh3 = (
        [x0 + i*dx for i in 0:nx-1], 
        [y0 + i*dy for i in 0:ny-1],
        [z0 + i*dz for i in 0:nz-1]
    )
    levelset3 = (x) -> (x[1] - a)^2 + (x[2] - a)^2 + (x[3] - a)^2 - R^2  # 3D sphere

    V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh3, levelset3)
    
    # Check that the results are the same
    @test length(V) == length(V1)
    @test length(cell_types) == length(cell_types1)
    @test length(interface_area) == length(Γ1)
    @test length(As) == length(A1)
    @test length(Ws) == length(W1)
    @test length(Bs) == length(B1)

    @test all(@. isapprox(V, V1, atol=1e-7))
    @test cell_types1 == cell_types
    @test all(@. isapprox(interface_area, Γ1, atol=1e-3))

    # Check face apertures (compare individual components)
    for i in 1:3
        @test all(@. isapprox(As[i], A1[i], atol=1e-5))
    end
    
    # Check staggered volumes
    for i in 1:3
        @test all(@. isapprox(Ws[i], W1[i], atol=1e-5))
    end
    
    # Check boundary fractions
    for i in 1:3
        @test all(@. isapprox(Bs[i], B1[i], atol=1e-5))
    end
    
    # Visualize results (optional)
    println("Volume comparison (first few cells):")
    println("Original: ", V[1:5])
    println("Implicit: ", V1[1:5])
    
    println("\nNumber of each cell type:")
    println("Empty cells: ", count(==(0), cell_types1))
    println("Full cells: ", count(==(1), cell_types1))
    println("Interface cells: ", count(==(-1), cell_types1))
    
    println("\nVolume fraction: ", sum(V1) / (nx*ny*nz*dx*dy*dz))
    println("Analytical volume fraction: ", 4π*R^3/(3*nx*ny*nz*dx*dy*dz))
end