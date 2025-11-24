using Test
using StaticArrays
using CartesianGeometry
using ImplicitIntegration
using LinearAlgebra

@testset "4D Geometry Integration" begin
    # Test with a simple 4D sphere centered at the origin
    # Φ(x,y,z,t) = x^2 + y^2 + z^2 + t^2 - r^2 (negative inside the sphere)
    r = 0.5  # radius
    Φ = x -> sum(x.^2) - r^2
    
    # Create a simple Cartesian grid
    n = 10
    mesh = (
        range(-1.0, 1.0, n+1) |> collect,  # x coordinates
        range(-1.0, 1.0, n+1) |> collect,  # y coordinates
        range(-1.0, 1.0, n+1) |> collect,  # z coordinates
        range(-1.0, 1.0, n+1) |> collect   # t coordinates
    )
    
    # Compute the geometric quantities
    V, cell_types, C_ω, C_γ, Γ, W, A, B = CartesianGeometry.implicit_integration(mesh, Φ)
    
    # Expected volume of a 4D sphere
    expected_volume = π^2 * r^4 / 2
    # Get the total computed volume
    computed_volume = sum(V)
    
    # Test if the computed volume is close to the expected volume
    @test isapprox(computed_volume, expected_volume, rtol=0.1)
    
    # Test surface area of a 4D sphere (3-sphere)
    # Surface area of a 3-sphere = 2π^2 * r^3
    expected_surface_area = 2 * π^2 * r^3
    computed_surface_area = sum(Γ)
    
    # Test if the computed surface area is close to the expected value
    @test isapprox(computed_surface_area, expected_surface_area, rtol=0.1)
    
    # Test cell types: should have some empty (0), full (1), and cut cells (-1)
    @test any(cell_types .== 0)  # some cells are empty
    @test any(cell_types .== 1)  # some cells are full
    @test any(cell_types .== -1) # some cells are cut
    
    
    # Analytical solution for a simple test case: translate sphere center slightly
    # This creates a well-defined inside and outside region
    offset = SVector(0.3, 0.3, 0.3, 0.3)
    Φ_off = x -> sum((x .- offset).^2) - r^2
    
    # Compute with offset sphere
    V1, cell_types1, C_ω1, C_γ1, Γ1, W1, A1, B1 = CartesianGeometry.implicit_integration(mesh, Φ_off)
    
    # Test staggered volumes W = (Wx, Wy, Wz, Wt)
    Ws = W
    # Check dimensions of W components
    @test length(W) == 4
    @test size(W[1]) == size(V)
    @test size(W[2]) == size(V)
    @test size(W[3]) == size(V)
    @test size(W[4]) == size(V)
    
    # Test face capacities A = (Ax, Ay, Az, At)
    As = A
    # Check dimensions of A components
    @test length(A) == 4
    @test size(A[1]) == size(V)
    @test size(A[2]) == size(V)
    @test size(A[3]) == size(V)
    @test size(A[4]) == size(V)
    
    # Test boundary fractions B = (Bx, By, Bz, Bt)
    Bs = B
    # Check dimensions of B components
    @test length(B) == 4
    @test size(B[1]) == size(V)
    @test size(B[2]) == size(V)
    @test size(B[3]) == size(V)
    @test size(B[4]) == size(V)
    
    # Test specific geometric relations:
    # 1. For full cells (cell_types = 1), W should equal dx*dy*dz*dt
    # 2. For empty cells (cell_types = 0), W should be 0
    dx = mesh[1][2] - mesh[1][1]
    dy = mesh[2][2] - mesh[2][1]
    dz = mesh[3][2] - mesh[3][1]
    dt = mesh[4][2] - mesh[4][1]
    cell_vol = dx * dy * dz * dt
    
    for i in eachindex(cell_types)
        if cell_types[i] == 1  # full cell
            @test isapprox(V[i], cell_vol, atol=1e-8)
        elseif cell_types[i] == 0  # empty cell
            @test isapprox(V[i], 0.0, atol=1e-8)
        end
    end

    # Test symmetry properties for a centered sphere
    # For a sphere centered at origin, volume should be symmetric in all directions
    center_idx = div.(size.(mesh), 2)
    
    # Get linear index of center cell and adjacent cells
    nx = length(mesh[1])-1
    ny = length(mesh[2])-1
    nz = length(mesh[3])-1
    nt = length(mesh[4])-1
    
    # Linear indexing function from the implementation
    linear_idx = (i, j, k, l) -> i + (j-1)*(nx+1) + (k-1)*(nx+1)*(ny+1) + (l-1)*(nx+1)*(ny+1)*(nz+1)
    
    # Check center and adjacent cells
    cx, cy, cz, ct = center_idx
    center_cell = linear_idx(cx, cy, cz, ct)
    pos_x = linear_idx(cx+1, cy, cz, ct)
    neg_x = linear_idx(cx-1, cy, cz, ct)
    pos_y = linear_idx(cx, cy+1, cz, ct)
    neg_y = linear_idx(cx, cy-1, cz, ct)
    pos_z = linear_idx(cx, cy, cz+1, ct)
    neg_z = linear_idx(cx, cy, cz-1, ct)
    pos_t = linear_idx(cx, cy, cz, ct+1)
    neg_t = linear_idx(cx, cy, cz, ct-1)
    
    # For a centered sphere, cells equidistant from center should have similar volumes
    if 1 <= pos_x <= length(V) && 1 <= neg_x <= length(V)
        @test isapprox(V[pos_x], V[neg_x], atol=1e-8)
    end
    
    if 1 <= pos_y <= length(V) && 1 <= neg_y <= length(V)
        @test isapprox(V[pos_y], V[neg_y], atol=1e-8)
    end
    
    if 1 <= pos_z <= length(V) && 1 <= neg_z <= length(V)
        @test isapprox(V[pos_z], V[neg_z], atol=1e-8)
    end
    
    if 1 <= pos_t <= length(V) && 1 <= neg_t <= length(V)
        @test isapprox(V[pos_t], V[neg_t], atol=1e-8)
    end
end