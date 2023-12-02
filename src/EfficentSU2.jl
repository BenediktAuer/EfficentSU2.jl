module EfficentSU2
    using StaticArrays
    import Base: *, similar, ones, adjoint,convert
    using LinearAlgebra
    import LinearAlgebra: tr, mul!,  adjoint!

    mutable struct SU2{T} <: FieldVector{2,T}
        z₁::T
        z₂::T
    end
    SU2(a,b) = SU2(@MVector [a,b])

    """
    getMatrix(a::SU2)

Returns the SU(2) Matrixrepresentation generatet by the two Elements of `a`

# Examples 
```julia-repl
julia> a = SU2(2+3im,1+4im)
2-element SU2{Complex{Int64}} with indices SOneTo(2):
 2 + 3im
 1 + 4im
julia> getMatrix(a)
2×2 StaticArraysCore.MMatrix{2, 2, Complex{Int64}, 4} with indices SOneTo(2)×SOneTo(2):
  2+3im  1+4im
 -1+4im  2-3im
```
"""
function getMatrix(a::SU2) 
        return @MMatrix [[a[1] a[2] ];[- conj(a[2]) conj(a[1])]]
    end
    # function Base.*(a::T,b::T) where T<:SU2
    #     return
    
    # end

*(a::SU2,b::SU2)= SU2(a[1]*b[1]-a[2]*conj(b[2]), a[1]*b[2]+a[2]*conj(b[1]))
tr(a::SU2) = a[1]+conj(a[1])

similar(a::SU2) = SU2(MArray{Tuple{2}}( Array{eltype(a)}(undef,2)))
"""
    mul!(res::SU2,a::SU2,b::SU2)
    Calculate a*b and store it in res
"""
function mul!(res::SU2,a::SU2,b::SU2)
    res[1] = a[1]*b[1]-a[2]*conj(b[2])
    res[2] = a[1]*b[2]+a[2]*conj(b[1])
    return
end
function ones(::Type{T}, dims::Tuple{Vararg{I, N}} where I<:Integer) where {T<:SU2,N}
    return reshape([ones(T) for i in 1:prod(dims)],dims)
end
"""
    renormalize!(a::SU2)
Normalizes the Matrix `a`such that  ``\\abs{z_1}^1+\\abs{z_2} =1``

```julia-repl
julia> a = SU2(2+3f0im,1+4f0im)
2-element SU2{ComplexF32} with indices SOneTo(2):
 2.0f0 + 3.0f0im
 1.0f0 + 4.0f0im
 julia> renormalize!(a)
 julia> a
 2-element SU2{ComplexF32} with indices SOneTo(2):
 0.36514837f0 + 0.5477225f0im
 0.18257418f0 + 0.73029673f0im
 julia> sum(abs2,a)
 0.9999999f0

```
"""
function renormalize!(a::SU2) 
     a./= norm(a)
    return
end
adjoint(a::SU2) = SU2(conj(a[1]),-a[2])
function adjoint!(res::SU2{T},a::SU2{T}) where T<:Number
        res[1] = conj(a[1])
        res[2] = -a[2]
        return
    end
    Base.promote_type(::Type{SU2{T}},::Type{SU2{S}}) where {T,S} = SU2{promote_type(T,S)}

    export SU2,getMatrix,*,tr,mul!,similar,ones,renormalize!,adjoint,adjoint!, promote_rule
end

