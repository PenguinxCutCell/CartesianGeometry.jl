using Test
using CartesianGeometry
using StaticArrays

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

    v_vofi, bary_vofi = CartesianGeometry.integrate(Tuple{0}, f, grid, Float64, zero; method=:vofi)
    v_jul, bary_jul = CartesianGeometry.integrate(Tuple{0}, f, grid, Float64, zero; method=:vofijul)

    @test isapprox(v_vofi[1], 0.5; atol=1e-6, rtol=1e-6)
    @test isapprox(v_jul[1], v_vofi[1]; atol=1e-6, rtol=1e-6)
    @test isapprox(bary_vofi[1][1], bary_jul[1][1]; atol=1e-4, rtol=1e-4)
    @test isapprox(bary_vofi[1][2], bary_jul[1][2]; atol=1e-4, rtol=1e-4)
end
