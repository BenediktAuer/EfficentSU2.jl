using EfficientSU2
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
    # mul!(c,a,b)
    # @test getMatrix(c) == @MMatrix([[-4+2im  4.0im];[4im -4-2im]])
end
@testset "Renorm and Adjoint" begin
    a = SU2(1+1f0*im,2f0*im)
   renormed =  renormalize(a)
    @test getMatrix(renormed*renormed')≈I #a' is the adjoint
end
# @testset "inplace multiplicattion" begin
#     a = SU2(1+1f0*im,2f0*im)
#     b = SU2(5+7f0*im,19f0*im)
#     c = SU2(-3+9f0*im,3f0*im)
#     res = similar(a)
#     resacc = similar(a)
#     mul!(res,resacc,a,b,c)
#     @test getMatrix(res) == getMatrix(a*b*c)
# end
# @testset "inplace add/sub" begin
#     a = SU2(1+1f0*im,2f0*im)
#     b = SU2(5+7f0*im,19f0*im)
#     res = similar(a)
#     add!(res,a,b)
#     @test getMatrix(res) == getMatrix(a+b)
#     sub!(res,a,b)
#     @test getMatrix(res) == getMatrix(a-b)
# end
#TODO add Test for ones check for deepcopys by modifing one entry and chekc if its the only one who has changed
