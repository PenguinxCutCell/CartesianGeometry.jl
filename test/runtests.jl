using CartesianGeometry
using Test
using StaticArrays

const T = Float64

@testset "Grid helper constructors" begin
    outer = 1:9
    inner = 1:8

    x = collocated(identity, outer, inner)
    @test first(x) == zero(eltype(x))
    @test last(x) == one(eltype(x))

    x = staggered(identity, outer, inner)
    @test 2 * length(inner) * first(x) == one(eltype(x))
end

@testset "1D hyper-sphere capacities" begin
    grid = collect(0.0:0.01:1.0)
    R = 0.2
    levelset = HyperSphere(R, (0.5,))

    V, bary, interface_length, cell_types, bary_interface = integrate(Tuple{0}, levelset, (grid,), T, nan)
    As = integrate(Tuple{1}, levelset, (grid,), T, nan)
    Ws = integrate(Tuple{0}, levelset, (grid,), T, nan, bary)
    Bs = integrate(Tuple{1}, levelset, (grid,), T, nan, bary)
    println(typeof(bary_interface))
    println(bary_interface)
    @test length(V) == length(bary) == length(interface_length) == length(cell_types) == length(bary_interface) == length(grid) 
    @test length(As) == length(Ws) == length(Bs)

    interface_lengths = filter(!isnan, interface_length)
    @test !isempty(interface_lengths)
    @test isapprox(sum(interface_lengths), 1.0; atol=1e-6)

    @test count(t -> t == -1, cell_types) > 0
    @test sum(V[.!isnan.(V)]) ≈ 2R atol=1e-6
end

@testset "2D hyper-sphere capacities" begin
    grid = (collect(0.0:0.05:1.0), collect(0.0:0.05:1.0))
    R = 0.4
    levelset = HyperSphere(R, (0.5, 0.5))

    V, bary, interface_length, cell_types, bary_interface = integrate(Tuple{0}, levelset, grid, T, nan)
    As = integrate(Tuple{1}, levelset, grid, T, nan)
    Ws = integrate(Tuple{0}, levelset, grid, T, nan, bary)
    Bs = integrate(Tuple{1}, levelset, grid, T, nan, bary)
    println(typeof(bary_interface))
    println(bary_interface)
    @test length(V) == length(bary) == length(interface_length) == length(cell_types) == length(bary_interface) == length(grid[1]) * length(grid[2])
    @test length(As) == length(Ws) == length(Bs) == length(grid)

    interface_lengths = filter(!isnan, interface_length)
    @test !isempty(interface_lengths)
    @test isapprox(sum(interface_lengths), 2 * π * R; atol=1e-6)

    @test count(t -> t == -1, cell_types) > 0
    @test sum(V[.!isnan.(V)]) ≈ π * R^2 atol=1e-6
end

@testset "3D hyper-sphere capacities" begin
    universe = (-1:11, -1:11, -1:11)
    node = (1:9, 1:9, 1:9)
    xyz = collocated.(identity, universe, node)
    R = 0.25
    levelset = HyperSphere(R, (0.5, 0.5, 0.5))

    V, bary, interface_area, cell_types, bary_interface = integrate(Tuple{0}, levelset, xyz, T, nan)
    As = integrate(Tuple{1}, levelset, xyz, T, nan)
    Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)
    Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)
    println(typeof(bary_interface))
    println(bary_interface)
    @test length(V) == length(bary) == length(interface_area) == length(cell_types) == length(bary_interface)
    @test length(As) == length(Ws) == length(Bs) == length(xyz)

    interface_area_vals = filter(!isnan, interface_area)
    @test !isempty(interface_area_vals)
    @test isapprox(sum(interface_area_vals), 4 * π * R ^ 2; atol=1e-3)

    @test typeof(Bs) == Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    @test count(t -> t == -1, cell_types) > 0
    @test sum(V[.!isnan.(V)]) ≈ (4 / 3) * π * R^3 atol=1e-4
end

@testset "4D hyper-sphere capacities" begin
    universe = (-1:11, -1:11, -1:11, -1:11)
    node = (1:9, 1:9, 1:9, 1:9)
    xyzw = collocated.(identity, universe, node)
    R = 0.5
    levelset = HyperSphere(R, (0.5, 0.5, 0.5, 0.5))

    V, bary, interface_hyperarea, cell_types, bary_interface = integrate(Tuple{0}, levelset, xyzw, T, nan; method=:vofijul)
    As = integrate(Tuple{1}, levelset, xyzw, T, nan; method=:vofijul)
    Ws = integrate(Tuple{0}, levelset, xyzw, T, nan, bary; method=:vofijul)
    Bs = integrate(Tuple{1}, levelset, xyzw, T, nan, bary; method=:vofijul)
    println(typeof(bary_interface))
    println(bary_interface)
    @test length(V) == length(bary) == length(interface_hyperarea) == length(cell_types) == length(bary_interface)
    @test length(As) == length(Ws) == length(Bs) == length(xyzw)

    interface_hyperarea_vals = filter(!isnan, interface_hyperarea)
    @test !isempty(interface_hyperarea_vals)
    @test isapprox(sum(interface_hyperarea_vals), 2 * π^2 * R ^ 3; atol=1e-1)

    @test typeof(Bs) == Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    @test count(t -> t == -1, cell_types) > 0
    @test sum(V[.!isnan.(V)]) ≈ (1 / 2) * π^2 * R^4 atol=1e-3

end


@testset "Vofi backends agree" begin
    f(x, y,_=0) = x + y - 1.0
    x = SVector(0.0, 1.0)
    y = SVector(0.0, 1.0)

    for backend in (:vofi, :vofijul)
        xex = zeros(Cdouble, 4)
        vol = CartesianGeometry.vofinit_dispatch!(backend, xex, f, x, y; nex=Cint.((1, 1)))

        @test isapprox(vol, 0.5; atol=1e-6, rtol=1e-6)
        @test isapprox(xex[1], 1 / 3; atol=1e-4, rtol=1e-4)
        @test isapprox(xex[2], 1 / 3; atol=1e-4, rtol=1e-4)
        @test CartesianGeometry.get_cell_type(f, x, y; backend=backend) == -1
    end
end

@testset "VofiJul 4D integration" begin
    f(x, y, z, w) = x + y + z + w - 1.0
    x = SVector(0.0, 1.0)
    y = SVector(0.0, 1.0)
    z = SVector(0.0, 1.0)
    w = SVector(0.0, 1.0)

    xex = zeros(Cdouble, 5)
    vol = CartesianGeometry.vofinit_dispatch!(:vofijul, xex, f, x, y, z, w; nex=Cint.((1, 1)))

    @test isapprox(vol, 1 / 24; atol=1e-6, rtol=1e-6)
    @test CartesianGeometry.get_cell_type(f, x, y, z, w; backend=:vofijul) == -1
end

@testset "integrate backend parity" begin
    f(x, y,_=0) = x + y - 1.0
    grid = (0.0:1.0:1.0, 0.0:1.0:1.0)

    v_vofi, bary_vofi, _, _, _ = CartesianGeometry.integrate(Tuple{0}, f, grid, Float64, zero; method=:vofi)
    v_jul, bary_jul, _, _, _ = CartesianGeometry.integrate(Tuple{0}, f, grid, Float64, zero; method=:vofijul)

    @test isapprox(v_vofi[1], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(v_jul[1], v_vofi[1]; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_vofi[1][1], bary_jul[1][1]; atol=1e-4, rtol=1e-4)
    @test isapprox(bary_vofi[1][2], bary_jul[1][2]; atol=1e-4, rtol=1e-4)
end

@testset "Interface centroid (linear level sets)" begin
    # 1D: interface at x = 0.5 inside the single cell
    grid1 = (0.0:1.0:1.0,)
    f1(x) = x - 0.5
    _, _, _, _, bary_if_1d = CartesianGeometry.integrate(Tuple{0}, f1, grid1, Float64, nan; method=:vofijul)
    @test isapprox(bary_if_1d[1][1], 0.5; atol=1e-6, rtol=1e-6)
    println(bary_if_1d)
    # 2D: plane x + y = 1 intersects square -> segment centroid at (0.5, 0.5)
    grid2 = (0.0:1.0:1.0, 0.0:1.0:1.0)
    f2(x, y) = x + y - 1.0
    _, _, _, _, bary_if_2d = CartesianGeometry.integrate(Tuple{0}, f2, grid2, Float64, nan; method=:vofijul)
    @test isapprox(bary_if_2d[1][1], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_2d[1][2], 0.5; atol=1e-6, rtol=1e-6)
    
    # 3D: plane x + y + z = 1 cuts cube -> triangular interface centroid at (1/3, 1/3, 1/3)
    grid3 = (0.0:1.0:1.0, 0.0:1.0:1.0, 0.0:1.0:1.0)
    f3(x, y, z) = x + y + z - 1.0
    _, _, _, _, bary_if_3d = CartesianGeometry.integrate(Tuple{0}, f3, grid3, Float64, nan; method=:vofijul)
    @test isapprox(bary_if_3d[1][1], 1 / 3; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_3d[1][2], 1 / 3; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_3d[1][3], 1 / 3; atol=1e-6, rtol=1e-6)

    # 4D: hyperplane x + y + z + w = 1 in hypercube -> simplex centroid at (1/4, 1/4, 1/4, 1/4)
    grid4 = (0.0:1.0:1.0, 0.0:1.0:1.0, 0.0:1.0:1.0, 0.0:1.0:1.0)
    f4(x, y, z, w) = x + y + z + w - 1.0
    _, _, _, _, bary_if_4d = CartesianGeometry.integrate(Tuple{0}, f4, grid4, Float64, nan; method=:vofijul)
    @test isapprox(bary_if_4d[1][1], 0.25; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_4d[1][2], 0.25; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_4d[1][3], 0.25; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_4d[1][4], 0.25; atol=1e-6, rtol=1e-6)

    # 4D Hypersphere: centroid at (0.5, 0.5, 0.5, 0.5)
    grid_sphere = (0.0:1.0:1.0, 0.0:1.0:1.0, 0.0:1.0:1.0, 0.0:1.0:1.0)
    f_sphere(x, y, z, w) = (x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2 + (w - 0.5)^2 - 0.25
    _, _, _, _, bary_if_sphere = CartesianGeometry.integrate(Tuple{0}, f_sphere, grid_sphere, Float64, nan; method=:vofijul)
    @test isapprox(bary_if_sphere[1][1], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_sphere[1][2], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_sphere[1][3], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_if_sphere[1][4], 0.5; atol=1e-6, rtol=1e-6)
end
