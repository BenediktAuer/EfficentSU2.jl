using EfficentSU2
using Test, StaticArrays,LinearAlgebra

@testset "Instantiation and getFields" begin
    a=SU2(1+1f0*im,2f0*im)
    @test typeof(a) <:FieldVector{2,ComplexF32}
    @test a[1] == 1+1f0im
    @test a[2] == 2f0im
end

@testset "getMatrix" begin
    a = SU2(2+3f0im,1+4f0im)
    @test getMatrix(a) == @MMatrix([[2.0+3.0im  1.0+4.0im];[-1+4im 2.0-3im]])
end
@testset "similar" begin
    a = SU2(2+3f0im,1+4f0im)
    @test typeof(a) == typeof(similar(a))
end

@testset "Multiplication" begin
    a = SU2(1+1f0*im,2f0*im)
    b = SU2(1+1f0*im,2f0*im)
    c = similar(a)
    @test getMatrix(a*b) == @MMatrix([[-4+2im  4.0im];[4im -4-2im]])
    mul!(c,a,b)
    @test getMatrix(c) == @MMatrix([[-4+2im  4.0im];[4im -4-2im]])
end

