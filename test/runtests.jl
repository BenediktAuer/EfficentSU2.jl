using EfficentSU2
using Test, StaticArrays

@testset "Instantiation and getFields" begin
    a=SU2(1+1f0*im,2f0*im)
    @test typeof(a) <:FieldVector{2,ComplexF32}
    @test a[1] == 1+1f0im
    @test a[2] == 2f0im
end
#TODO Add test for getMatrix()
#TODO Add test for Base.*