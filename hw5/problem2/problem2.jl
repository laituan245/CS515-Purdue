using LinearAlgebra
using DelimitedFiles
using SparseArrays

# ========================= Helper functions =========================
function extract_diagonal(A)
    # Compute D and D_inverse from A
    index = 1
    I, J, V, V_inverse = zeros(A.n), zeros(A.n), zeros(A.n), zeros(A.n)
    for j=1:length(A.colptr)-1 # for each column ...
        for nzi=A.colptr[j]:A.colptr[j+1]-1 # for each entry in the column
            i = A.rowval[nzi]
            v = A.nzval[nzi]
            if i == j
                I[index] = i
                J[index] = j
                V[index] = v
                V_inverse[index] = 1.0 / v
                index += 1
            end
        end
    end

    D = sparse(I, J, V, A.n, A.n)
    D_inverse = sparse(I, J, V_inverse, A.n, A.n)
    return D, D_inverse
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
function jacobi_method(P, A, b)
    # A is a CSC matrix
    # b is a vector
    A, b = P * A, P * b

    D, D_inverse = extract_diagonal(A)
    N = A - D

    c = D_inverse * b
    h = D_inverse * N

    iterations = 0
    x = rand(A.n)
    while norm(csc_matvec(A, x) - b) / norm(b) > 1e-4
        x = c - csc_matvec(h, x)
        iterations += 1
    end
    return x, iterations
end
jacobi_x, jacobi_iterations = jacobi_method(I, T'-I, -b)
println(jacobi_x[140])
println(jacobi_x[47])
println(jacobi_iterations)
