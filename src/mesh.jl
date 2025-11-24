"""
    collocated(f, outer, inner)

Created a one-dimensional mesh where the first and last
elements (excluding halo) are collocated with the boundary.

"""
function collocated(f, outer, inner)
    x = Vector{Float64}(undef, length(outer))

    lo = findfirst(isequal(first(inner)), outer)
    up = findfirst(isequal(last(inner)), outer)

    for i in eachindex(outer)
        x[i] = f((i-lo) / (up-lo+1))
    end

    x
#    CVector(x, tuple(outer))
end

"""
    staggered(f, outer, inner)

Created a one-dimensional mesh where the first and last
elements (excluding halo) are staggered with the boundary.

"""
function staggered(f, outer, inner)
    x = Vector{Float64}(undef, length(outer))

    lo = findfirst(isequal(first(inner)), outer)
    up = findfirst(isequal(last(inner)), outer)

    for i in eachindex(outer)
        x[i] = f((2(i-lo) + 1) / 2(up-lo+1))
    end

    x
#    CVector(x, tuple(outer))
end
