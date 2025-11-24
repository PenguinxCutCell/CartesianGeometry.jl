function vofinit!(xex, f, x::Number; nex=Cint.((1, 1)), backend=:vofi)
    t = f(x)

    val = one(x)

    isnonpositive(t) && return val
    isnonnegative(t) && return zero(val)

    nothing
end

"""

1d implementation complies with Vofi 2.0's API (see 2d and 3d).

!!! note "Convention"

    Even if cell is empty, return full centroid coordinates.

"""
function vofinit!(xex, f, x::SVector; nex=Cint.((1, 1)), backend=:vofi)
    t = SVector{2}(f(i) for i in x)

    val = x[2] - x[1]

    if all(isnonpositive, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        isone(last(nex)) && (xex[end] = zero(xex[1]))
        return val
    end

    if all(isnonnegative, t)
        isone(first(nex)) && (xex[1] = sum(x) / 2)
        isone(last(nex)) && (xex[end] = zero(xex[1]))
        return zero(val)
    end

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    if isnonnegative(t[1])
        isone(first(nex)) && (xex[1] = (ξ + x[2]) / 2)
        isone(last(nex)) && (xex[end] = 1.0)
        return x[2] - ξ
    end

    if isnonnegative(t[2])
        isone(first(nex)) && (xex[1] = (x[1] + ξ) / 2)
        isone(last(nex)) && (xex[end] = 1.0)
        return ξ - x[1]
    end

    nothing
end

function vofinit!(xex, f, x::SVector{2}, y::Number; nex=Cint.((1, 1)), backend=:vofi)
    t = SVector{2}(f(i, y) for i in x)

    val = x[2] - x[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (x[2] * t[1] - x[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return x[2] - ξ
    isnonnegative(t[2]) && return ξ - x[1]

    nothing
end

function vofinit!(xex, f, x::Number, y::SVector{2}; nex=Cint.((1, 1)), backend=:vofi)
    t = SVector{2}(f(x, j) for j in y)

    val = y[2] - y[1]

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    ξ = (y[2] * t[1] - y[1] * t[2]) / (t[1] - t[2])

    isnonnegative(t[1]) && return y[2] - ξ
    isnonnegative(t[2]) && return ξ - y[1]

    nothing
end

"""

Call Vofi 2.0 for exact integration. ###

"""
function vofinit!(xex, f, x::SVector, y::SVector; nex=Cint.((1, 1)), backend=:vofi)
    t = SMatrix{2,2}(f(i, j) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    val * vofi_getcc(backend, f, x0, h0, xex, 2; nex=nex)
end

function vofinit!(xex, f, x::Number, y::SVector{2}, z::SVector{2}; nex=Cint.((1, 1)), backend=:vofi)
    t = SMatrix{2,2}(f(x, j, k) for j in y, k in z)

    val = (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], zero(val)))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], one(val)))

    nex = Cint.((0, 0))
    val * vofi_getcc(backend, (yval, zval) -> f(x, yval, zval), x0, h0, xex, 2; nex=nex)
end

function vofinit!(xex, f, x::SVector{2}, y::Number, z::SVector{2}; nex=Cint.((1, 1)), backend=:vofi)
    t = SMatrix{2,2}(f(i, y, k) for i in x, k in z)

    val = (z[2] - z[1]) * (x[2] - x[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((z[1], x[1], zero(val)))
    h0 = Cdouble.((z[2]-z[1], x[2]-x[1], one(val)))

    nex = Cint.((0, 0))
    val * vofi_getcc(backend, (zval, xval) -> f(xval, y, zval), x0, h0, xex, 2; nex=nex)
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::Number; nex=Cint.((1, 1)), backend=:vofi)
    t = SMatrix{2,2}(f(i, j, z) for i in x, j in y)

    val = (x[2] - x[1]) * (y[2] - y[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], zero(val)))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], one(val)))

    nex = Cint.((0, 0))
    val * vofi_getcc(backend, (xval, yval) -> f(xval, yval, z), x0, h0, xex, 2; nex=nex)
end

"""

Call Vofi 2.0 for exact integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector, z::SVector; nex=Cint.((1, 1)), backend=:vofi)
    t = SArray{Tuple{2,2,2}}(f(i, j, k) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    val * vofi_getcc(backend, f, x0, h0, xex, 3; nex=nex)
end

"""

Call Vofi extended for exact 4D integration.

"""
function vofinit!(xex, f, x::SVector, y::SVector, z::SVector, w::SVector; nex=Cint.((1, 1)), backend=:vofi)
    t = SArray{Tuple{2,2,2,2}}(f(i, j, k, l) for i in x, j in y, k in z, l in w)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1]) * (w[2] - w[1])

    if all(isnonpositive, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[4] = sum(w) / 2
        xex[end] = zero(xex[1])
        return val
    end

    if all(isnonnegative, t)
        xex[1] = sum(x) / 2
        xex[2] = sum(y) / 2
        xex[3] = sum(z) / 2
        xex[4] = sum(w) / 2
        xex[end] = zero(xex[1])
        return zero(val)
    end

    x0 = Cdouble.((x[1], y[1], z[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1], w[2]-w[1]))

    val * vofi_getcc(backend, f, x0, h0, xex, 4; nex=nex)
end

# Fonctions pour les combinaisons d'arguments 4D avec une dimension fixe
function vofinit!(xex, f, x::Number, y::SVector{2}, z::SVector{2}, w::SVector{2}; nex=Cint.((0, 0)), backend=:vofi)
    t = SArray{Tuple{2,2,2}}(f(x, j, k, l) for j in y, k in z, l in w)

    val = (y[2] - y[1]) * (z[2] - z[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((y[1], z[1], w[1]))
    h0 = Cdouble.((y[2]-y[1], z[2]-z[1], w[2]-w[1]))

    val * vofi_getcc(backend, (yval, zval, wval) -> f(x, yval, zval, wval), x0, h0, xex, 3; nex=nex)
end

function vofinit!(xex, f, x::SVector{2}, y::Number, z::SVector{2}, w::SVector{2}; nex=Cint.((0, 0)), backend=:vofi)
    t = SArray{Tuple{2,2,2}}(f(i, y, k, l) for i in x, k in z, l in w)

    val = (x[2] - x[1]) * (z[2] - z[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], z[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], z[2]-z[1], w[2]-w[1]))

    val * vofi_getcc(backend, (xval, zval, wval) -> f(xval, y, zval, wval), x0, h0, xex, 3; nex=nex)
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::Number, w::SVector{2}; nex=Cint.((0, 0)), backend=:vofi)
    t = SArray{Tuple{2,2,2}}(f(i, j, z, l) for i in x, j in y, l in w)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (w[2] - w[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], w[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], w[2]-w[1]))

    val * vofi_getcc(backend, (xval, yval, wval) -> f(xval, yval, z, wval), x0, h0, xex, 3; nex=nex)
end

function vofinit!(xex, f, x::SVector{2}, y::SVector{2}, z::SVector{2}, w::Number; nex=Cint.((0, 0)), backend=:vofi)
    t = SArray{Tuple{2,2,2}}(f(i, j, k, w) for i in x, j in y, k in z)

    val = (x[2] - x[1]) * (y[2] - y[1]) * (z[2] - z[1])

    all(isnonpositive, t) && return val
    all(isnonnegative, t) && return zero(val)

    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2]-x[1], y[2]-y[1], z[2]-z[1]))

    val * vofi_getcc(backend, (xval, yval, zval) -> f(xval, yval, zval, w), x0, h0, xex, 3; nex=nex)
end

const DEFAULT_NPT = Cint.((4, 4, 4, 4))
const DEFAULT_NVIS = Cint.((0, 0))

function normalize_vofi_backend(method)
    if method === :vofijul || method === VofiJul
        return :vofijul
    elseif method === :vofi || method === Vofinit || method === nothing
        return :vofi
    else
        throw(ArgumentError("Unknown VOF backend: $method"))
    end
end

function _call_level_set(f, coords, ndim)
    args = ntuple(i -> coords[i], ndim)
    if applicable(f, args...)
        return f(args...)
    elseif applicable(f, coords...)
        return f(coords...)
    elseif applicable(f, coords)
        return f(coords)
    else
        throw(ArgumentError("Level set function is not callable with $ndim dimensions"))
    end
end

_wrap_vofijul_integrand(f, ndim) = coords -> _call_level_set(f, coords, ndim)

function vofi_getcc(backend, f, x0, h0, xex, ndim; nex=Cint.((1, 1)),
                    npt=DEFAULT_NPT, nvis=DEFAULT_NVIS)
    if backend === :vofijul
        return VofiJul.vofi_get_cc(_wrap_vofijul_integrand(f, ndim), nothing,
                                   x0, h0, xex, Int.(nex), Int.(npt), Int.(nvis), Int(ndim))
    elseif backend === :vofi
        return Vofinit.getcc(f, x0, h0, xex, Cint(ndim);
                             nex=Cint.(nex), npt=Cint.(npt), nvis=Cint.(nvis))
    end
    throw(ArgumentError("Unsupported VOF backend: $backend"))
end

function vofi_get_cell_type(backend, f, x0, h0, ndim)
    if backend === :vofijul
        return VofiJul.vofi_get_cell_type(_wrap_vofijul_integrand(f, ndim),
                                          nothing, x0, h0, Int(ndim))
    elseif backend === :vofi
        return Vofinit.getcelltype(f, x0, h0, Cint(ndim))
    end
    throw(ArgumentError("Unsupported VOF backend: $backend"))
end

compute_cell_type(backend, f, args...) = get_cell_type(f, args...; backend=backend)

function get_cell_type(f, x::SVector{2}; backend=:vofi)
    if backend === :vofijul
        x0 = Cdouble.((x[1], 0.0, 0.0))
        h0 = Cdouble.((x[2] - x[1], 1.0, 1.0))
        return vofi_get_cell_type(backend, f, x0, h0, 1)
    elseif backend === :vofi
        t = SVector{2}(f(i) for i in x)
        if all(isnonpositive, t)
            return 1.0
        end
        if all(isnonnegative, t)
            return 0.0
        end
        if isnonnegative(t[1])
            return -1.0
        end
        if isnonnegative(t[2])
            return -1.0
        end
        return nothing
    end
    throw(ArgumentError("Unsupported VOF backend for cell type: $backend"))
end

function get_cell_type(f, x::SVector{2}, y::SVector{2}; backend=:vofi)
    x0 = Cdouble.((x[1], y[1], 0.0))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], 1.0))
    return vofi_get_cell_type(backend, f, x0, h0, 2)
end

function get_cell_type(f, x::SVector{2}, y::SVector{2}, z::SVector{2}; backend=:vofi)
    x0 = Cdouble.((x[1], y[1], z[1]))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], z[2] - z[1]))
    return vofi_get_cell_type(backend, f, x0, h0, 3)
end

function get_cell_type(f, x::SVector{2}, y::SVector{2}, z::SVector{2}, w::SVector{2}; backend=:vofi)
    x0 = Cdouble.((x[1], y[1], z[1], w[1]))
    h0 = Cdouble.((x[2] - x[1], y[2] - y[1], z[2] - z[1], w[2] - w[1]))
    return vofi_get_cell_type(backend, f, x0, h0, 4)
end



# Ajouter une fonction de dispatch à la fin du fichier:
"""
    vofinit_dispatch!(method, xex, f, args...; nex=Cint.((1, 1)))

Dispatch helper selecting either Vofi C (`:vofi`) or VofiJul (`:vofijul`).
"""
function vofinit_dispatch!(method, xex, f, args...; nex=Cint.((1, 1)))
    backend = normalize_vofi_backend(method)
    return vofinit!(xex, f, args...; nex=nex, backend=backend)
end
