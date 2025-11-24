struct Dirichlet{F,T,R}
    f::F
    xyz::T
    ranges::R

    Dirichlet(f::F, xyz::T, ranges::R) where {F,T,R} =
        new{F,T,R}(f, xyz, ranges)
    Dirichlet(f::Type{F}, xyz::T, ranges::R) where {F,T,R} =
        new{Type{F},T,R}(f, xyz, ranges)
end

function apply!(v, border::Dirichlet)
    (; f, xyz, ranges) = border
    (; outer, inner) = ranges

    rv = reshape(v, length.(outer)...)

    rxyz = map(xyz, outer) do x, r
        reshape(x, length(r))
    end

    inner = findin.(inner, outer)
    outer = eachindex.(outer)

    _apply!(rv, f, rxyz, outer, inner)

    v
end

@generated function _apply!(v::ArrayAbstract{N}, f, xyz, outer, inner) where {N}
    quote
        @ntuple($N, x) = xyz

        for index in EdgeIterator(outer, inner)
            @ntuple($N, i) = Tuple(index)
            @nref($N, v, i) = @ncall($N, f, d -> x_d[i_d])
        end

        v
    end
end

#=
struct Dirichlet2{F}
    f::F

    Dirichlet2(f::F) where {F} = new{F}(f)
    Dirichlet2(f::Type{F}) where {F} = new{Type{F}}(f)
end

homogeneous(args...) = zero(mapreduce(eltype, promote_type, args))
measure(args::SVector{2}...) = mapreduce(*, args) do el
    el[2] - el[1]
end

=#
