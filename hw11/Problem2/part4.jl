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
