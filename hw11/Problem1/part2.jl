using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("hw11_functions.jl")

# Load data
data = readdlm("poisson2D-data.csv",',')
A = sparse(Int.(data[:,1]), Int.(data[:,2]),(data[:,3]))
A = (A + A')./2
b = vec(readdlm("poisson2D-rhs.csv"))

# Define matrix C
sigma = 1.7 * 1e-2
C = A - sigma * I

# Use preconditioner + GMRES
Lchol = ichol(C)
M1 = LowerPrecond(ichol(C))
x_pgm, hist_pgm = gmres(C, b, restart = 30, tol = 1.0e-6 * norm(b), maxiter = 10000, Pl=M1, log=true)
println(hist_pgm)
rel_residuals_pgm = log10.(hist_pgm.data[:resnorm] / norm(b))
savefig(plot(rel_residuals_pgm, xlabel="iterations", ylabel="Log Relative residuals",
             label="preconditioner + GMRES", linewidth=2), "preconditioner_gmres.png")


# Use preconditioner + MINRES
using Krylov # ] clone https://github.com/JuliaSmoothOptimizers/Krylov.jl.git
using LinearOperators
P1 = opInverse(LowerTriangular(ichol(C)))
Mprecond = P1'*P1;
x_pmin, hist_pmin = Krylov.minres(C, b, M=Mprecond, itmax=10000, rtol=1e-10)
rel_residuals_pmin = log10.(hist_pmin.residuals / norm(b))
savefig(plot(rel_residuals_pmin, xlabel="iterations", ylabel="Log Relative residuals",
             label="preconditioner + MINRES", linewidth=2), "preconditioner_minres.png")
