using LinearAlgebra
using DelimitedFiles
using SparseArrays

# ========================= Helper functions =========================
function extract_diagonal(A)
    # Compute A_diagonal
    A_diagonal = zeros(A.m)
    for j=1:length(A.colptr)-1 # for each column ...
        for nzi=A.colptr[j]:A.colptr[j+1]-1 # for each entry in the column
            i = A.rowval[nzi]
            v = A.nzval[nzi]
            if i == j
                A_diagonal[i] = v
            end
        end
    end
    return A_diagonal
end

# Input the transition probability matrix T
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

# Prepare the vector b
b = ones((140, 1))
b[134] = 0

# Solve for correct x
x =  (T' - I) \ -b
