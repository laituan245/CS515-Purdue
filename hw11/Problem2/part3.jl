using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

# Implement a function for solving the Poisson's equation using the Jacobi method
function solve_poisson_using_jacobi(n)
    F = poisson_setup(n,n, (x,y) -> 1)
    X = zeros(size(F)...)
    while norm(poisson_residual(X,F)) > 1e-10
        X = apply_poisson_jacobi(X, F)
    end
    return X
end

X_poisson = solve_poisson_using_jacobi(31)
X_direct = solve_poisson_direct(poisson_setup(31,31, (x,y) -> 1))
println(norm(X_poisson - X_direct))

# Implement code to analyze the error trend
function analyze_error_trend(n, iters=1000)
    errors = zeros(iters)
    F = poisson_setup(n,n, (x,y) -> 1)
    X_correct = solve_poisson_direct(F)

    X_jacobi = zeros(size(F)...)
    for i=1:iters
        X_jacobi = apply_poisson_jacobi(X_jacobi, F)
        errors[i] = norm(X_jacobi - X_correct) / norm(X_correct)
    end

    ratios = zeros(iters)
    for i=2:iters
        ratios[i] = errors[i] / errors[i-1]
    end
    savefig(plot(2:iters, ratios[2:iters], linewidth=2, label="Error Ratio"), "jacobi_error_trend.png")
    println(string("last ratio = ", ratios[iters]))
end
analyze_error_trend(31)
