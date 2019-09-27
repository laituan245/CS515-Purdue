using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Input the transition probability matrix T
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

# Prepare the vector b
b = ones((140, 1))
b[134] = 0

# Solve for correct x
x =  (T' - I) \ -b
