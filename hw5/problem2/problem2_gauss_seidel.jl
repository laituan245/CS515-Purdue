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

function extract_L_and_U(A)
    # Compute L (the strict lower part) and U (the strict upper part) of A
    U_I, U_J, U_V = Int64[], Int64[], Float64[]
    L_I, L_J, L_V = Int64[], Int64[], Float64[]
    for j=1:length(A.colptr)-1 # for each column ...
        for nzi=A.colptr[j]:A.colptr[j+1]-1 # for each entry in the column
            i = A.rowval[nzi]
            v = A.nzval[nzi]
            if j > i
                push!(U_I, i)
                push!(U_J, j)
                push!(U_V, v)
            end
            if i > j
                push!(L_I, i)
                push!(L_J, j)
                push!(L_V, v)
            end
        end
    end

    L = sparse(L_I, L_J, L_V, A.n, A.n)
    U = sparse(U_I, U_J, U_V, A.n, A.n)
    return L, U
end

function sparse_transpose(A)
    # Compute A^T
    I, J, V = Int64[], Int64[], Float64[]
    for j=1:length(A.colptr)-1 # for each column ...
        for nzi=A.colptr[j]:A.colptr[j+1]-1 # for each entry in the column
            i = A.rowval[nzi]
            v = A.nzval[nzi]
            push!(I, j)
            push!(J, i)
            push!(V, v)
        end
    end
    return sparse(I, J, V, A.n, A.n)
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

# ========================= Gauss-Seidel Implementation =========================
function gauss_seidel_method(P, A, b)
    # A is a CSC matrix
    # b is a vector
    A, b = P * A, P * b

    D, D_inverse = extract_diagonal(A)
    L, U = extract_L_and_U(A)

    c1 = D_inverse * b
    c2 = D_inverse * U
    c3 = sparse_transpose(D_inverse * L)

    iterations = 0
    x = rand(A.n)
    while norm(csc_matvec(A, x) - b) / norm(b) > 1e-4
        h_k = c1 - csc_matvec(c2, x)
        for i in 1:A.n
            x[i] = h_k[i]
            for nzi=c3.colptr[i]:c3.colptr[i+1]-1 # for each entry in the column
                x[i] -= x[c3.rowval[nzi]] * c3.nzval[nzi]
            end
        end
        iterations += 1
    end
    return x, iterations
end

gauss_seidel_x, gauss_seidel_iterations = gauss_seidel_method(I, T'-I, -b)
println(gauss_seidel_x[140])
println(gauss_seidel_x[47])
println(gauss_seidel_iterations)
