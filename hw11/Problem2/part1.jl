using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

function solve(n)
    println(string("n = ", n))
    @time solve_poisson_direct(poisson_setup(n,n, (x,y) -> 1))
    println("End\n")
end

solve(31)
solve(63)
solve(127)
solve(255)
solve(511)
solve(1023)
