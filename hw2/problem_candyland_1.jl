using LinearAlgebra
using DelimitedFiles
using SparseArrays

# ---------------------------------- PART 2 ----------------------------------
# # Input the transition probability matrix T
# data = readdlm("candyland-matrix.csv",',')
# TI = Int.(data[:,1])
# TJ = Int.(data[:,2])
# TV = data[:,3]
# T = sparse(TI,TJ, TV, 140,140)
#
# # Prepare the vector b
# b = ones((140, 1))
# b[134] = 0
#
# # Solve for x
# x =  (T' - I) \ -b
# println(x[140])
# println(x[47])
# ----------------------------------------------------------------------------

# ---------------------------------- PART 3 ----------------------------------
# Input the transition probability matrix T
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

function calculate_expected_length(T)
    p = T[:, 140]
    expected_length, k = 0, 0
    while norm(p) > 1e-18
        k += 1
        expected_length += k * p[134]
        p = T * p
    end
    return expected_length
end
println(calculate_expected_length(T))
# ----------------------------------------------------------------------------
