"""
!!! warning

    [Vofinit](https://github.com/vlc1/Vofinit.jl) requires
    level set function to take 3 positional arguments.

"""
struct HyperSphere{N,T}
    radius::T
    center::NTuple{N,T}
end

function (object::HyperSphere{1})(x, _...)
    (; radius, center) = object
    (x - center[1]) ^ 2 - radius ^ 2
end

function (object::HyperSphere{2})(x, y, _...)
    (; radius, center) = object
    (x - center[1]) ^ 2 + (y - center[2]) ^ 2 - radius ^ 2
end

function (object::HyperSphere{3})(x, y, z)
    (; radius, center) = object
    (x - center[1]) ^ 2 + (y - center[2]) ^ 2 + (z - center[3]) ^ 2 - radius ^ 2
end

function (object::HyperSphere{4})(x, y, z, w)
    (; radius, center) = object
    (x - center[1]) ^ 2 + (y - center[2]) ^ 2 + (z - center[3]) ^ 2 + (w - center[4]) ^ 2 - radius ^ 2
end
