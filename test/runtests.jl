using CartesianGeometry
using Test

@testset "CartesianGeometry.jl" begin
    outer = 1:9
    inner = 1:8

    x = collocated(identity, outer, inner)
    @test first(x) == zero(eltype(x))
    @test last(x) == one(eltype(x))

    x = staggered(identity, outer, inner)
    @test 2length(inner) * first(x) == one(eltype(x))
end

include("vofi_backends.jl")
