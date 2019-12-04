using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Distributions
include("hw11_functions.jl")
pyplot()

# Load data
data = readdlm("poisson2D-data.csv",',')
A = sparse(Int.(data[:,1]), Int.(data[:,2]),(data[:,3]))
A = (A + A')./2
b = vec(readdlm("poisson2D-rhs.csv"))

# Define matrix C
sigma = 1.7 * 1e-2
C = A - sigma * I

# Show the eigenvalues of the matrix before preconditioning
before_matrix = Array(C)
before_eigvals = eigvals(before_matrix)
savefig(histogram(before_eigvals, bins=100, xlabel="Eigenvalues", ylabel="Counts"), "before_preconditioning.png")

# Show the eigenvalues of the matrix after preconditioning
using Krylov # ] clone https://github.com/JuliaSmoothOptimizers/Krylov.jl.git
using LinearOperators
M1 = LowerPrecond(ichol(C))
after_matrix = Array(M1.lower * C)
after_eigvals = eigvals(after_matrix)
savefig(histogram(after_eigvals, bins=100, xlabel="Eigenvalues", ylabel="Counts"), "after_preconditioning.png")
