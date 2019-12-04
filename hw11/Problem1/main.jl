using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers

# Load data
data = readdlm("poisson2D-data.csv",',')
A = sparse(Int.(data[:,1]), Int.(data[:,2]),(data[:,3]))
A = (A + A')./2
b = vec(readdlm("poisson2D-rhs.csv"))

# Do MINRES and GMRES
sigma = 1.7 * 1e-2
x_min, hist_min = minres(A - sigma * I, b, tol = 1e-6, maxiter = 10000, log=true)
x_gm, hist_gm = gmres(A - sigma * I, b, restart = 30, maxiter = 10000, tol = 1e-6, log=true)

# Calculate the relative residuals
rel_residuals_min = hist_min.data[:resnorm] / norm(b)
rel_residuals_gm = hist_gm.data[:resnorm] / norm(b)
