#=

          ?           ?           ?           ?           ?


y_4 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  V_{1,3}  |  V_{2,3}  |  V_{3,3}  |  V_{4,3}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_3 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  V_{1,2}  |  V_{2,2}  |  V_{3,2}  |  V_{4,2}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_2 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  V_{1,1}  |  V_{2,1}  |  V_{3,1}  |  V_{4,1}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_1 +-----------+-----------+-----------+-----------+
   x_1         x_2         x_3         x_4         x_5

              (and same for barycenters...)

=#

"""

    integrate(Tuple{0}, f, xyz, T, bc)

Computes volume-specific (`Tuple{0}`) apertures of the first kind.

Returns a `Tuple` where

1. The first component is a `Vector{T}` that stores the wet volumes of each cell,
1. The second component is a `Vector{SVector{N,T}}` that stores the coordinates of the wet barycenters.

Wherever the moments can not be computed, the function `bc` is applied to the element type.

# Arguments

- `f`: the level set function.
- `xyz`: the Cartesian coordinates of the lattice nodes.
- `T`: the `eltype` of the moments.
- `bc`: the boundary condition (*e.g.* `nan` or `zero`).

!!! note

    To simplify the computation of second-kind moments, barycenters of `-f` are stored in empty cells.

"""
function integrate(::Type{Tuple{0}}, f, xyz, T, bc; method=:vofi)
    N, len = length(xyz), prod(length.(xyz))

    v = Vector{T}(undef, len)
    bary = Vector{SVector{N,T}}(undef, len)
    interface_norm = Vector{T}(undef, len)
    cell_types = Vector{T}(undef, len) 

    integrate!((v, bary, interface_norm, cell_types), Tuple{0}, f, xyz, bc; method=method)
end

# ND volume
@generated function integrate!(moms, ::Type{Tuple{0}},
                              f, xyz::NTuple{N}, bc; method=:vofi) where {N}
    quote
        # Mêmes axes et indices...
        input = only.(axes.(xyz))
        output = droplast.(input)
        linear = LinearIndices(input)
        cartesian = CartesianIndices(output)

        (v, bary, interface_norm, cell_types) = moms
        @ntuple($N, x) = xyz
        xex = zeros(Cdouble, 4)
        backend = normalize_vofi_backend(method)

        for index in cartesian
            n = linear[index]
            @ntuple($N, i) = Tuple(index)
            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))
            
            # Utiliser la méthode spécifiée
            v[n] = @ncall($N, vofinit_dispatch!, backend, xex, f, y)
            bary[n] = @ncall($N, SVector, d -> xex[d])
            interface_norm[n] = xex[end]
            
            # Déterminer le type de cellule
            cell_types[n] = @ncall($N, compute_cell_type, backend, f, y)
        end

        # boundary conditions
        halo = EdgeIterator(CartesianIndices(input), cartesian)

        for index in halo
            n = linear[index]

            v[n] = bc(eltype(v))
            bary[n] = bc(eltype(bary))
            interface_norm[n] = bc(eltype(interface_norm))
            cell_types[n] = bc(eltype(cell_types))
        end

        return moms
    end
end

#=

    ?           ?           ?           ?           ?


y_4 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
 A_{1,3}     A_{2,3}     A_{3,3}     A_{4,3}     A_{5,3}
    |           |           |           |           |
    |           |           |           |           |
y_3 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
 A_{1,2}     A_{2,2}     A_{3,2}     A_{4,2}     A_{5,2}
    |           |           |           |           |
    |           |           |           |           |
y_2 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
 A_{1,1}     A_{2,1}     A_{3,1}     A_{4,1}     A_{5,1}
    |           |           |           |           |
    |           |           |           |           |
y_1 +-----------+-----------+-----------+-----------+
   x_1         x_2         x_3         x_4         x_5

=#
#=

y_4 +--A_{1,4}--+--A_{2,4}--+--A_{3,4}--+--A_{4,4}--+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_3 +--A_{1,3}--+--A_{2,3}--+--A_{3,3}--+--A_{4,3}--+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_2 +--A_{1,2}--+--A_{2,2}--+--A_{3,2}--+--A_{4,2}--+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_1 +--A_{1,1}--+--A_{2,1}--+--A_{3,1}--+--A_{4,1}--+     ?
   x_1         x_2         x_3         x_4         x_5


=#

"""

    integrate(Tuple{1}, f, xyz, T, bc)

Computes area-specific (`Tuple{1}`) apertures of the first kind.

Returns a `NTuple` where each element corresponds to  direction.

# Arguments

- `f`: the level set function.
- `xyz`: the Cartesian coordinates of the lattice nodes.
- `T`: the `eltype` of the moments.
- `bc`: the boundary condition (*e.g.* `nan` or `zero`).

"""
function integrate(::Type{Tuple{1}}, f, xyz, T, bc; method=:vofi)
    len = prod(length.(xyz))

    moms = map(xyz) do _
        Vector{T}(undef, len)
    end

    integrate!(moms, Tuple{1}, f, xyz, bc; method=method)
end

# 1D surface
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{1}, bc; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    for n in linear
        moms[1][n] = vofinit_dispatch!(method, xex, f, x[n])
    end

    return moms
end

# 2D surface
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{2}, bc; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = input[1], droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        moms[1][n] = vofinit_dispatch!(method, xex, f, x[i],
                                      SVector(y[j], y[j+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), input[2]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        moms[2][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      y[j])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    return moms
end

# 3D surface
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{3}, bc; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = input[1], droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        moms[1][n] = vofinit_dispatch!(method, xex, f, x[i],
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), input[2], droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        moms[2][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      y[j],
                                      SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    # z faces

    output = droplast(input[1]), droplast(input[2]), input[3]
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        moms[3][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      z[k])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[3][n] = bc(eltype(moms[3]))
    end

    return moms
end
