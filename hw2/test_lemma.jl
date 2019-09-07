using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Input the transition probability matrix T
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

A = Matrix(T)
S = (I-A)^(-2)
println(ndims(S))

S2 = I
for i in 2:10000
    global S2
    S2 = S2 + i * A^(i-1)
end
println("test 1")
println(S[10, 1:10])
println(S2[10, 1:10])
println("\ntest 2")
println(S[20, 10:20])
println(S2[20, 10:20])
