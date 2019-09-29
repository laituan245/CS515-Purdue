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

function csc_matvec(A, x)
    # Compute the matrix-vector product Ax
      y = zeros(A.n)
      for j=1:length(A.colptr)-1 # for each column ...
          for nzi=A.colptr[j]:A.colptr[j+1]-1 # for each entry in the column
              i = A.rowval[nzi]
              v = A.nzval[nzi]
              y[i] += v*x[j]
          end
      end
      return y
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

# ========================= JACOBI Implementation =========================
function jacobi_method(A, b)
    # A is a CSC matrix
    # b is a vector
    iterations = 0
    x = rand(A.n)
    A_diagonals = extract_diagonal(A)
    while norm(A*x - b) / norm(b) > 1e-4
        x = (b - A*x + x .* A_diagonals) ./ A_diagonals
        iterations += 1
    end
    return x, iterations
end
jacobi_x, jacobi_iterations = jacobi_method(T'-I, -b)
