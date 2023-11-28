module EfficentSU2
    using StaticArrays
    struct SU2{T<:Number}
        m::MVector{2, T}
    end
    SU2(a,b) = SU2(@MVector [a,b])
    export SU2
end
