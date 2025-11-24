"""

integrate(T::Type{<:Tuple}, f, xyz, S, bc, bary)

Computes volume- (`T=Tuple{0}`) and surface-specific (`T=Tuple{1}`) apertures of the second kind.

# Arguments

- `f`: Level set function.
- `xyz::NTuple{N}`: the Cartesian coordinates of the lattice nodes.
- `bc`: Boundary condition.
- `bary`: Barycenters.

!!! warning

    The last barycenters along each dimension are undefined.

"""
function integrate(T::Type{<:Tuple}, f, xyz, S, bc, bary; method=:vofi)

    moms = map(xyz) do _
        similar(bary, S)
    end

    integrate!(moms, T, f, xyz, bc, bary; method=method)
end

#=

    ?           ?           ?           ?           ?


y_4 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    ?        W_{2,3}     W_{3,3}     W_{4,3}        ?
    |           |           |           |           |
    |           |           |           |           |
y_3 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    ?        W_{2,2}     W_{3,2}     W_{4,2}        ?
    |           |           |           |           |
    |           |           |           |           |
y_2 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    ?        W_{2,1}     W_{3,1}     W_{4,1}        |
    |           |           |           |           |
    |           |           |           |           |
y_1 +-----------+-----------+-----------+-----------+
   x_1         x_2         x_3         x_4         x_5

=#
#=

y_4 +-----?-----+-----?-----+-----?-----+-----?-----+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_3 +--W_{1,3}--+--W_{2,3}--+--W_{3,3}--+--W_{4,3}--+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_2 +--W_{1,2}--+--W_{2,2}--+--W_{3,2}--+--W_{4,2}--+     ?
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
    |           |           |           |           |
y_1 +-----?-----+-----?-----+-----?-----+-----?-----+     ?
   x_1         x_2         x_3         x_4         x_5


=#

# 1D volume
function integrate!(moms, ::Type{Tuple{0}},
                    f, xyz::NTuple{1}, bc, bary; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = (dropends(input[1]),)
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, = Tuple(index)

        left, right = i-1, n

        moms[1][n] = vofinit_dispatch!(method, xex, f,
                              SVector(bary[left][1], bary[right][1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    return moms
end

# 2D volume
function integrate!(moms, ::Type{Tuple{0}},
                    f, xyz::NTuple{2}, bc, bary; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = dropends(input[1]), droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        left, right = linear[i-1, j], n

        moms[1][n] = vofinit_dispatch!(method, xex, f,
                              SVector(bary[left][1], bary[right][1]),
                              SVector(y[j], y[j+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), dropends(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        left, right = linear[i, j-1], n

        moms[2][n] = vofinit_dispatch!(method, xex, f,
                              SVector(x[i], x[i+1]),
                              SVector(bary[left][2], bary[right][2]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    return moms
end

# 2D volume
function integrate!(moms, ::Type{Tuple{0}},
                    f, xyz::NTuple{3}, bc, bary; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = dropends(input[1]), droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        left, right = linear[i-1, j, k], n

        moms[1][n] = vofinit_dispatch!(method, xex, f,
                              SVector(bary[left][1], bary[right][1]),
                              SVector(y[j], y[j+1]),
                              SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), dropends(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        left, right = linear[i, j-1, k], n

        moms[2][n] = vofinit_dispatch!(method, xex, f,
                              SVector(x[i], x[i+1]),
                              SVector(bary[left][2], bary[right][2]),
                              SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    # z faces

    output = droplast(input[1]), droplast(input[2]), dropends(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        left, right = linear[i, j, k-1], n

        moms[3][n] = vofinit_dispatch!(method, xex, f,
                              SVector(x[i], x[i+1]),
                              SVector(y[j], y[j+1]),
                              SVector(bary[left][3], bary[right][3]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[3][n] = bc(eltype(moms[3]))
    end

    return moms
end


#=

          ?           ?           ?           ?           ?


y_4 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  B_{1,3}  |  B_{2,3}  |  B_{3,3}  |  B_{4,3}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_3 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  B_{1,2}  |  B_{2,2}  |  B_{3,2}  |  B_{4,2}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_2 +-----------+-----------+-----------+-----------+
    |           |           |           |           |
    |           |           |           |           |
    |  B_{1,1}  |  B_{2,1}  |  B_{3,1}  |  B_{4,1}  |     ?
    |           |           |           |           |
    |           |           |           |           |
y_1 +-----------+-----------+-----------+-----------+
   x_1         x_2         x_3         x_4         x_5

                   (same for x and y)

=#

# 1D surface
function integrate!(moms, ::Type{Tuple{1}},
                    f, xyz::NTuple{1}, bc, bary; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = (droplast(input[1]),)
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, = Tuple(index)

        moms[1][n] = vofinit_dispatch!(method, xex, f, bary[i][1])
    end

    return moms
end

# 2D surface
function integrate!(moms, ::Type{Tuple{1}},
                    f, xyz::NTuple{2}, bc, bary; method=:vofi)

    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = droplast(input[1]), droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        moms[1][n] = vofinit_dispatch!(method, xex, f,
                              bary[n][1],
                              SVector(y[j], y[j+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), droplast(input[2])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j = Tuple(index)

        moms[2][n] = vofinit!(xex, f,
                              SVector(x[i], x[i+1]),
                              bary[n][2])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    return moms
end

# 3D surface
function integrate!(moms, ::Type{Tuple{1}},
                    f, xyz::NTuple{3}, bc, bary; method=:vofi)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = droplast(input[1]), droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)
        moms[1][n] = vofinit_dispatch!(method, xex, f,
                              bary[n][1],
                              SVector(y[j], y[j+1]),
                              SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces

    output = droplast(input[1]), droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        moms[2][n] = vofinit_dispatch!(method, xex, f,
                              SVector(x[i], x[i+1]),
                              bary[n][2],
                              SVector(z[k], z[k+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    # z faces

    output = droplast(input[1]), droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k = Tuple(index)

        moms[3][n] = vofinit_dispatch!(method, xex, f,
                              SVector(x[i], x[i+1]),
                              SVector(y[j], y[j+1]),
                              bary[n][3])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[3][n] = bc(eltype(moms[3]))
    end

    return moms
end