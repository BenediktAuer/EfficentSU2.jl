module EfficentSU2
    using StaticArrays
    import Base.*,Base.similar, LinearAlgebra.norm,LinearAlgebra.tr, LinearAlgebra.mul!
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
 """
    *(a::SU2,b::SU2)
#TODO Add DFocstring
TBW
"""
*(a::SU2,b::SU2)= SU2(a[1]*b[1]-a[2]*conj(b[2]), a[1]*b[2]+a[2]*conj(b[1]))
tr(a::SU2) = a[1]+conj(a[1])

similar(a::SU2) = SU2(MArray{Tuple{2}}( Array{eltype(a)}(undef,2)))
function mul!(res::SU2,a::SU2,b::SU2)
    res[1] = a[1]*b[1]-a[2]*conj(b[2])
    res[2] = a[1]*b[2]+a[2]*conj(b[1])
    return
end
    export SU2,getMatrix,*,tr,mul!,similar
end

