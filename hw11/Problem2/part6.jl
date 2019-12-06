using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

# Directly solving a Poisson problem where n = 31
X_31_true = solve_poisson_direct(poisson_setup(31,31, (x,y) -> 1))

# Directly solving a Poisson problem where n = 63
X_63_true = solve_poisson_direct(poisson_setup(63,63, (x,y) -> 1))

# Interpolate the solution from n = 31 to a n = 63 solution
X_63_interpolate = interpolate(X_31_true)
