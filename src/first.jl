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
    bary_interface = Vector{SVector{N,T}}(undef, len)

    integrate!((v, bary, interface_norm, cell_types, bary_interface), Tuple{0}, f, xyz, bc; method=method)
end

"""
    integrate_centroid(f, xyz, T, bc; method=:vofi)

Compute the interface centroid in each cell (when available from the backend) by
leveraging `integrate(Tuple{0}, ...)`. Falls back to `bc` wherever the centroid
cannot be computed (e.g. with the C backend).
"""
function integrate_centroid(f, xyz, T, bc; method=:vofi)
    _, _, _, _, bary_interface = integrate(Tuple{0}, f, xyz, T, bc; method=method)
    return bary_interface
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

        (v, bary, interface_norm, cell_types, bary_interface) = moms
        @ntuple($N, x) = xyz
        xex = zeros(Cdouble, 4)
        backend = normalize_vofi_backend(method)

        for index in cartesian
            n = linear[index]
            @ntuple($N, i) = Tuple(index)
            @nextract($N, y, d -> SVector(x_d[i_d], x_d[i_d+1]))
            coords = @ntuple($N, y)
            
            # Utiliser la méthode spécifiée
            v[n] = @ncall($N, vofinit_dispatch!, backend, xex, f, y)
            bary[n] = @ncall($N, SVector, d -> xex[d])
            interface_norm[n] = xex[end]
            
            # Déterminer le type de cellule
            cell_types[n] = @ncall($N, compute_cell_type, backend, f, y)

            centroid = interface_centroid(backend, f, coords...)
            if centroid === nothing
                bary_interface[n] = bc(eltype(bary_interface))
            else
                bary_interface[n] = centroid
            end
        end

        # boundary conditions
        halo = EdgeIterator(CartesianIndices(input), cartesian)

        for index in halo
            n = linear[index]

            v[n] = bc(eltype(v))
            bary[n] = bc(eltype(bary))
            interface_norm[n] = bc(eltype(interface_norm))
            cell_types[n] = bc(eltype(cell_types))
            bary_interface[n] = bc(eltype(bary_interface))
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

# 4D surface
function integrate!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{4}, bc; method=:vofijul)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z, w = xyz
    xex = zeros(Cdouble, 4)

    # x faces

    output = input[1], droplast(input[2]), droplast(input[3]), droplast(input[4])
    cartesian = CartesianIndices(output)

    for index in cartesian
        n = linear[index]
        i, j, k, l = Tuple(index)

        moms[1][n] = vofinit_dispatch!(method, xex, f, x[i],
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]),
                                      SVector(w[l], w[l+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]

        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces
    output = droplast(input[1]), input[2], droplast(input[3]), droplast(input[4])
    cartesian = CartesianIndices(output)
    for index in cartesian
        n = linear[index]
        i, j, k, l = Tuple(index)

        moms[2][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      y[j],
                                      SVector(z[k], z[k+1]),
                                      SVector(w[l], w[l+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]

        moms[2][n] = bc(eltype(moms[2]))
    end

    # z faces
    output = droplast(input[1]), droplast(input[2]), input[3], droplast(input[4])
    cartesian = CartesianIndices(output)
    for index in cartesian
        n = linear[index]
        i, j, k, l = Tuple(index)
        moms[3][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      z[k],
                                      SVector(w[l], w[l+1]))
    end
    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[3][n] = bc(eltype(moms[3]))
    end
    
    # w faces
    output = droplast(input[1]), droplast(input[2]), droplast(input[3]), input[4]
    cartesian = CartesianIndices(output)
    for index in cartesian
        n = linear[index]
        i, j, k, l = Tuple(index)
        moms[4][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]),
                                      w[l])
    end
    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[4][n] = bc(eltype(moms[4]))
    end

    return moms
end

#####################################################################
# THREADED VERSIONS
#####################################################################

"""
    integrate_threaded(Tuple{0}, f, xyz, T, bc; method=:vofi)

Multithreaded version of `integrate(Tuple{0}, ...)`.
Computes volume-specific (`Tuple{0}`) apertures of the first kind using multiple threads.

See [`integrate`](@ref) for more details.
"""
function integrate_threaded(::Type{Tuple{0}}, f, xyz, T, bc; method=:vofi)
    N, len = length(xyz), prod(length.(xyz))

    v = Vector{T}(undef, len)
    bary = Vector{SVector{N,T}}(undef, len)
    interface_norm = Vector{T}(undef, len)
    cell_types = Vector{T}(undef, len)
    bary_interface = Vector{SVector{N,T}}(undef, len)

    integrate_threaded!((v, bary, interface_norm, cell_types, bary_interface), Tuple{0}, f, xyz, bc; method=method)
end

# 1D volume - threaded version
function integrate_threaded!(moms, ::Type{Tuple{0}},
                              f, xyz::NTuple{1}, bc; method=:vofi)
    input = only.(axes.(xyz))
    output = droplast.(input)
    linear = LinearIndices(input)
    cartesian = CartesianIndices(output)

    (v, bary, interface_norm, cell_types, bary_interface) = moms
    x, = xyz
    backend = normalize_vofi_backend(method)

    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, = Tuple(index)
        y1 = SVector(x[i], x[i+1])
        coords = (y1,)
        
        xex = zeros(Cdouble, 4)
        
        v[n] = vofinit_dispatch!(backend, xex, f, y1)
        bary[n] = SVector(xex[1])
        interface_norm[n] = xex[end]
        
        cell_types[n] = compute_cell_type(backend, f, y1)

        centroid = interface_centroid(backend, f, coords...)
        if centroid === nothing
            bary_interface[n] = bc(eltype(bary_interface))
        else
            bary_interface[n] = centroid
        end
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]
        v[n] = bc(eltype(v))
        bary[n] = bc(eltype(bary))
        interface_norm[n] = bc(eltype(interface_norm))
        cell_types[n] = bc(eltype(cell_types))
        bary_interface[n] = bc(eltype(bary_interface))
    end

    return moms
end

# 2D volume - threaded version
function integrate_threaded!(moms, ::Type{Tuple{0}},
                              f, xyz::NTuple{2}, bc; method=:vofi)
    input = only.(axes.(xyz))
    output = droplast.(input)
    linear = LinearIndices(input)
    cartesian = CartesianIndices(output)

    (v, bary, interface_norm, cell_types, bary_interface) = moms
    x1, x2 = xyz
    backend = normalize_vofi_backend(method)

    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i1, i2 = Tuple(index)
        y1 = SVector(x1[i1], x1[i1+1])
        y2 = SVector(x2[i2], x2[i2+1])
        coords = (y1, y2)
        
        xex = zeros(Cdouble, 4)
        
        v[n] = vofinit_dispatch!(backend, xex, f, y1, y2)
        bary[n] = SVector(xex[1], xex[2])
        interface_norm[n] = xex[end]
        
        cell_types[n] = compute_cell_type(backend, f, y1, y2)

        centroid = interface_centroid(backend, f, coords...)
        if centroid === nothing
            bary_interface[n] = bc(eltype(bary_interface))
        else
            bary_interface[n] = centroid
        end
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]
        v[n] = bc(eltype(v))
        bary[n] = bc(eltype(bary))
        interface_norm[n] = bc(eltype(interface_norm))
        cell_types[n] = bc(eltype(cell_types))
        bary_interface[n] = bc(eltype(bary_interface))
    end

    return moms
end

# 3D volume - threaded version
function integrate_threaded!(moms, ::Type{Tuple{0}},
                              f, xyz::NTuple{3}, bc; method=:vofi)
    input = only.(axes.(xyz))
    output = droplast.(input)
    linear = LinearIndices(input)
    cartesian = CartesianIndices(output)

    (v, bary, interface_norm, cell_types, bary_interface) = moms
    x1, x2, x3 = xyz
    backend = normalize_vofi_backend(method)

    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i1, i2, i3 = Tuple(index)
        y1 = SVector(x1[i1], x1[i1+1])
        y2 = SVector(x2[i2], x2[i2+1])
        y3 = SVector(x3[i3], x3[i3+1])
        coords = (y1, y2, y3)
        
        xex = zeros(Cdouble, 4)
        
        v[n] = vofinit_dispatch!(backend, xex, f, y1, y2, y3)
        bary[n] = SVector(xex[1], xex[2], xex[3])
        interface_norm[n] = xex[end]
        
        cell_types[n] = compute_cell_type(backend, f, y1, y2, y3)

        centroid = interface_centroid(backend, f, coords...)
        if centroid === nothing
            bary_interface[n] = bc(eltype(bary_interface))
        else
            bary_interface[n] = centroid
        end
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]
        v[n] = bc(eltype(v))
        bary[n] = bc(eltype(bary))
        interface_norm[n] = bc(eltype(interface_norm))
        cell_types[n] = bc(eltype(cell_types))
        bary_interface[n] = bc(eltype(bary_interface))
    end

    return moms
end

# 4D volume - threaded version
function integrate_threaded!(moms, ::Type{Tuple{0}},
                              f, xyz::NTuple{4}, bc; method=:vofijul)
    input = only.(axes.(xyz))
    output = droplast.(input)
    linear = LinearIndices(input)
    cartesian = CartesianIndices(output)

    (v, bary, interface_norm, cell_types, bary_interface) = moms
    x1, x2, x3, x4 = xyz
    backend = normalize_vofi_backend(method)

    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i1, i2, i3, i4 = Tuple(index)
        y1 = SVector(x1[i1], x1[i1+1])
        y2 = SVector(x2[i2], x2[i2+1])
        y3 = SVector(x3[i3], x3[i3+1])
        y4 = SVector(x4[i4], x4[i4+1])
        coords = (y1, y2, y3, y4)
        
        xex = zeros(Cdouble, 5)
        
        v[n] = vofinit_dispatch!(backend, xex, f, y1, y2, y3, y4)
        bary[n] = SVector(xex[1], xex[2], xex[3], xex[4])
        interface_norm[n] = xex[end]
        
        cell_types[n] = compute_cell_type(backend, f, y1, y2, y3, y4)

        centroid = interface_centroid(backend, f, coords...)
        if centroid === nothing
            bary_interface[n] = bc(eltype(bary_interface))
        else
            bary_interface[n] = centroid
        end
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)

    for index in halo
        n = linear[index]
        v[n] = bc(eltype(v))
        bary[n] = bc(eltype(bary))
        interface_norm[n] = bc(eltype(interface_norm))
        cell_types[n] = bc(eltype(cell_types))
        bary_interface[n] = bc(eltype(bary_interface))
    end

    return moms
end

"""
    integrate_threaded(Tuple{1}, f, xyz, T, bc; method=:vofi)

Multithreaded version of `integrate(Tuple{1}, ...)`.
Computes area-specific (`Tuple{1}`) apertures of the first kind using multiple threads.

See [`integrate`](@ref) for more details.
"""
function integrate_threaded(::Type{Tuple{1}}, f, xyz, T, bc; method=:vofi)
    len = prod(length.(xyz))

    moms = map(xyz) do _
        Vector{T}(undef, len)
    end

    integrate_threaded!(moms, Tuple{1}, f, xyz, bc; method=method)
end

# 1D surface - threaded version
function integrate_threaded!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{1}, bc; method=:vofi)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, = xyz

    # Use @threads for parallel execution
    Threads.@threads for n in collect(linear)
        xex = zeros(Cdouble, 4)
        moms[1][n] = vofinit_dispatch!(method, xex, f, x[n])
    end

    return moms
end

# 2D surface - threaded version
function integrate_threaded!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{2}, bc; method=:vofi)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y = xyz

    # x faces
    output = input[1], droplast(input[2])
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j = Tuple(index)
        xex = zeros(Cdouble, 4)

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
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j = Tuple(index)
        xex = zeros(Cdouble, 4)

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

# 3D surface - threaded version
function integrate_threaded!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{3}, bc; method=:vofi)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z = xyz

    # x faces
    output = input[1], droplast(input[2]), droplast(input[3])
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k = Tuple(index)
        xex = zeros(Cdouble, 4)

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
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k = Tuple(index)
        xex = zeros(Cdouble, 4)

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
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k = Tuple(index)
        xex = zeros(Cdouble, 4)

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

# 4D surface - threaded version
function integrate_threaded!(moms, ::Type{Tuple{1}}, f, xyz::NTuple{4}, bc; method=:vofijul)
    input = only.(axes.(xyz))
    linear = LinearIndices(input)

    x, y, z, w = xyz

    # x faces
    output = input[1], droplast(input[2]), droplast(input[3]), droplast(input[4])
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k, l = Tuple(index)
        xex = zeros(Cdouble, 4)

        moms[1][n] = vofinit_dispatch!(method, xex, f, x[i],
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]),
                                      SVector(w[l], w[l+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[1][n] = bc(eltype(moms[1]))
    end

    # y faces
    output = droplast(input[1]), input[2], droplast(input[3]), droplast(input[4])
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k, l = Tuple(index)
        xex = zeros(Cdouble, 4)

        moms[2][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      y[j],
                                      SVector(z[k], z[k+1]),
                                      SVector(w[l], w[l+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[2][n] = bc(eltype(moms[2]))
    end

    # z faces
    output = droplast(input[1]), droplast(input[2]), input[3], droplast(input[4])
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k, l = Tuple(index)
        xex = zeros(Cdouble, 4)

        moms[3][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      z[k],
                                      SVector(w[l], w[l+1]))
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[3][n] = bc(eltype(moms[3]))
    end
    
    # w faces
    output = droplast(input[1]), droplast(input[2]), droplast(input[3]), input[4]
    cartesian = CartesianIndices(output)
    cartesian_arr = collect(cartesian)

    Threads.@threads for idx in eachindex(cartesian_arr)
        index = cartesian_arr[idx]
        n = linear[index]
        i, j, k, l = Tuple(index)
        xex = zeros(Cdouble, 4)

        moms[4][n] = vofinit_dispatch!(method, xex, f, SVector(x[i], x[i+1]),
                                      SVector(y[j], y[j+1]),
                                      SVector(z[k], z[k+1]),
                                      w[l])
    end

    halo = EdgeIterator(CartesianIndices(input), cartesian)
    for index in halo
        n = linear[index]
        moms[4][n] = bc(eltype(moms[4]))
    end

    return moms
end
