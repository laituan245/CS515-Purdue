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

# Evaluate the error
error = norm(X_63_interpolate - X_63_true)
println(string("The error is ", error))

# Make a mesh-plot
function make_mesh_plot(mat)
    m, n = size(mat)
    x = zeros(m * n); y = zeros(m * n); z = zeros(m * n);
    index = 1;
    for i=1:m
        for j=1:n
            x[index] = i
            y[index] = j
            z[index] = mat[i,j]
            index += 1
        end
    end
    return surface(x, y, z)
end

diff = abs.(X_63_interpolate - X_63_true)
pyplot()
savefig(make_mesh_plot(diff), "mesh_plot.png")
