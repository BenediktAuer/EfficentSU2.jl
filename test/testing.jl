using EfficientSU2,LinearAlgebra,Random
function get_randomSU2!(R,ϵ)
    #distance = 1/2--1/2 gives interval -1/2 to 1/2
    getRandUniformly!(R,1)
    r24 = view(R,2:4)
    norms = norm(r24)
       r24 .*= ϵ/norm(r24)
     R[1] = sign(R[1])*sqrt(1-ϵ^2)
     v = SU2(R[1]+R[4]*im,R[3]+R[2]*im)
     return v

end
"""
    getRandUniformly(::Type{T},distance,nums::Int64) where T
    returns `nums` RAndomNumbers Uniformly in the Intervall 0-distance to 0+distance 
TBW
"""
function getRandUniformly(::Type{T},distance,nums::Int64) where T
    rand(T,nums)./distance .-T(distance/2)
end

function getRandUniformly!(x::T,distance) where {N<:Number,T<:Array{N}}
    rand!(x)
    x.= x ./distance .-distance/2
end
a = ones(Float64,4)
b = get_randomSU2!(a,0.1)
c = get_randomSU2!(a,0.1)
get_randomSU2!(a,0.1)
get_randomSU2!(a,0.1)
get_randomSU2!(a,0.1)
@show c
