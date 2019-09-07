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
b_star = ones((140, 1))
b_star[5] = 0; b_star[135] = 0
b_star[35] = 0; b_star[136] = 0
b_star[47] = 0; b_star[86] = 0; b_star[117] = 0
b_star[134] = 0


# Solve for x
x =  (T' - I) \ -b_star
println(x[140])
