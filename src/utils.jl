const ArrayAbstract{N,T} = AbstractArray{T,N}

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonnegative(x) = x ≥ zero(x)
isnonpositive(x) = x ≤ zero(x)

droplast(this::Base.OneTo) = Base.OneTo(length(this)-1)
droplast(this::AbstractRange) = this[begin:end-1]

dropends(this::AbstractRange) = this[begin+1:end-1]

for T in [Float16, Float32, Float64]
    @eval nan(::Type{$T}) = $T(NaN)
    @eval nan(::Type{SVector{N,$T}}) where {N} =
        SVector{N,$T}(ntuple(Returns(nan($T)), N))
end
