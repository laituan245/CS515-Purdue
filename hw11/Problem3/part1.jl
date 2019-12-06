using DelimitedFiles
using LinearAlgebra

# Read the eigenfaces matrix and compute svd using Julia built-in function
A = readdlm(download("http://www.cs.purdue.edu/homes/dgleich/cs515-2019/homeworks/Yale_64.csv"),',')
M = A'
U,Sig,V = svd(M)
println(Sig)
