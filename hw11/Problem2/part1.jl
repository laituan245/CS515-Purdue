using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("multigrid_functions.jl")

function solve_poisson(n)
    return @timed solve_poisson_direct(poisson_setup(n,n, (x,y) -> 1))
end


measured_times = zeros(6)
measured_times[1] = log(2, solve_poisson(31)[2])
measured_times[2] = log(2, solve_poisson(63)[2])
measured_times[3] = log(2, solve_poisson(127)[2])
measured_times[4] = log(2, solve_poisson(255)[2])
measured_times[5] = log(2, solve_poisson(511)[2])
measured_times[6] = log(2, solve_poisson(1023)[2])

ns = zeros(6)
ns[1] = log(2, 31)
ns[2] = log(2, 63)
ns[3] = log(2, 127)
ns[4] = log(2, 255)
ns[5] = log(2, 511)
ns[6] = log(2, 1023)

savefig(plot(ns, measured_times, xlabel="Log2(n)", ylabel="Log2(Times measured in seconds)", linewidth=2, label="Time"), "time_analysis.png")
