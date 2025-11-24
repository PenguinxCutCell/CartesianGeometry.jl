using Test
using CartesianGeometry
using StaticArrays
using ImplicitIntegration

const T = Float64

@testset "2D Geometry Integration" begin
    # Setup mesh and level set
    R = 0.21
    a = 0.5
    nx, ny = 40, 40
    dx, dy = 1.0/nx, 1.0/ny
    x0, y0 = 0.0, 0.0
    mesh = ([x0 + i*dx for i in 0:nx], [y0 + i*dy for i in 0:ny])

    # Use a circle as the level set function
    levelset = (x, y, _=0) -> (x - a)^2 + (y - a)^2 - R^2

    # Compute moments using CartesianGeometry's integrate
    V, bary, interface_length, cell_types = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero)
    As = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero)
    Ws = CartesianGeometry.integrate(Tuple{0}, levelset, mesh, T, zero, bary)
    Bs = CartesianGeometry.integrate(Tuple{1}, levelset, mesh, T, zero, bary)

    # Now test the implicit integration implementation
    levelset2 = (x) -> (x[1] - a)^2 + (x[2] - a)^2 - R^2  # 2D circle as implicit function

    V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = implicit_integration(mesh, levelset2)
    
    # Check dimensions match
    @test length(V) == length(V1)
    @test length(cell_types) == length(cell_types1)
    @test length(interface_length) == length(Γ1)
    
    # Check tuple sizes for vectors
    @test length(As) == length(A1)
    @test length(Ws) == length(W1)
    @test length(Bs) == length(B1)

    # Check values match (with tolerance)
    @test all(@. isapprox(V, V1, atol=1e-7))
    @test cell_types1 == cell_types
    @test all(@. isapprox(interface_length, Γ1, atol=1e-5))

    # Check face apertures (compare individual components)
    for i in 1:2
        @test all(@. isapprox(As[i], A1[i], atol=1e-1))
    end
    
    # Check staggered volumes
    for i in 1:2
        @test all(@. isapprox(Ws[i], W1[i], atol=1e-5))
    end
    
    # Check boundary fractions
    for i in 1:2
        @test all(@. isapprox(Bs[i], B1[i], atol=1e-1))
    end
    
    # Print some statistics for debugging
    println("Volume comparison (first few cells):")
    println("Original: ", V[1:5])
    println("Implicit: ", V1[1:5])
    
    println("\nNumber of each cell type:")
    println("Empty cells: ", count(==(0), cell_types1))
    println("Full cells: ", count(==(1), cell_types1))
    println("Interface cells: ", count(==(-1), cell_types1))
    
    # Calculate area of the circle and compare with analytical value
    circle_area = sum(V1)
    analytical_area = π * R^2
    println("\nTotal area from integration: ", circle_area)
    println("Analytical circle area: ", analytical_area)
    println("Relative error: ", abs(circle_area - analytical_area)/analytical_area)
end