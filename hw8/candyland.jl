using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Input the transition probability matrix T
data = readdlm("data/candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

# Float64
function calculate_expected_length(T)
    p = T[:, 140]
    expected_length, k = 0.0, 0.0
    while norm(p) > 1e-18
        k += 1
        expected_length += k * p[134]
        p = T * p
    end
    return expected_length
end
println(calculate_expected_length(T))
# ----------------------------------------------------------------------------
