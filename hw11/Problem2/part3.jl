using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

function solve_poisson_using_jacobi(n)
    F = poisson_setup(n,n, (x,y) -> 1)
    X = zeros(size(F)...)
    while norm(poisson_residual(X,F)) > 1e-10
        X = apply_poisson_jacobi(X, F)
    end
    return X
end

X_poisson = solve_poisson_using_jacobi(63)
X_direct = solve_poisson_direct(poisson_setup(63,63, (x,y) -> 1))
println(norm(X_poisson - X_direct))
