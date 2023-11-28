using EfficentSU2
using Test, StaticArrays

@testset "EfficentSU2.jl" begin
    @test typeof(SU2(1*im,2*im).m) == MVector{2,Complex{Int64}}
end
