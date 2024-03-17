using EfficientSU2,LinearAlgebra,Random,BenchmarkTools   
a = SU2(1+1f0*im,2f0*im)
b = SU2(5+7f0*im,19f0*im)
c = SU2(-3+9f0*im,3f0*im)
res = similar(a)
resacc = similar(a)
@benchmark begin
    for _ in 1:2_560_000

    mul!($res,$resacc,$a,$b,$c) 
    end
end

@benchmark $a*$b*$c evals = 2_560_000


getMatrix(res) == getMatrix(a*b*c)

# VSCodeServer.@profview begin 
#     for _ in 1:100_000_000
#     mul!(res,resacc,a,b,c)
#     end
# end
@benchmark $a*$b*$c
@benchmark     mul!(res,resacc,$a,$b,$c) 
@benchmark $a.+$b
a.+b
a+b
