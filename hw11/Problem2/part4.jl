using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

function apply_poisson_gauss_seidel(X,B)
  Y = zeros(size(X)...)
  for j=2:size(B,2)-1
    for i=2:size(B,1)-1
      Y[i,j] = B[i,j]/4 +(X[i+1,j] + Y[i-1,j] + Y[i,j-1] + X[i,j+1])/4
    end
  end
  return Y
end

# Implement code to analyze the error trend
function analyze_error_trend(n, iters=1500)
    errors = zeros(iters)
    F = poisson_setup(n,n, (x,y) -> 1)
    X_correct = solve_poisson_direct(F)

    X_gauss_seidel = zeros(size(F)...)
    for i=1:iters
        X_gauss_seidel = apply_poisson_gauss_seidel(X_gauss_seidel, F)
        errors[i] = norm(X_gauss_seidel - X_correct) / norm(X_correct)
    end

    ratios = zeros(iters)
    for i=2:iters
        ratios[i] = errors[i] / errors[i-1]
    end
    savefig(plot(2:iters, ratios[2:iters], linewidth=2, label="Error Ratio (Gauss Seidel)"), "gauss_seidel_error_trend.png")
    println(string("last ratio = ", ratios[iters]))
    println(string("last error = ", errors[iters]))

end
analyze_error_trend(31)
