using LinearAlgebra
using DelimitedFiles
using SparseArrays

# ------------- Coordinate Descent / Gradient Descent Algorithms -------------
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

function update_g(old_g, A, A_diagonal, index)
  # old_g is the gradient vector from previous iteration
  # A is a CSC matrix
  # A_diagonal is a vector containing the values of the diagonal of A
  # index is the selected coordinate
  new_g = old_g
  colptr = A.colptr; rowval = A.rowval; nzval = A.nzval
  scalar = old_g[index] / A_diagonal[index]
  for nzi=colptr[index]:colptr[index+1]-1
      new_g[rowval[nzi]] -= (scalar * nzval[nzi])
  end
  return new_g
end

function cyclic_coordinate_descent(A, b)
    total_work = 0
    A_diagonal = extract_diagonal(A)
    x = rand(A.n)
    g = A * x - b
    i = 1
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-10
        x[i] = x[i] - (g[i] / A_diagonal[i])
        g = update_g(g, A, A_diagonal, i)
        total_work += (A.colptr[i+1]-A.colptr[i])
        relative_residual = norm(g) / norm(b)
        if (i == A.n) i = 1 else i += 1 end
    end
    return x, total_work
end

function random_coordinate_descent(A, b)
    total_work = 0
    A_diagonal = extract_diagonal(A)
    x = rand(A.n)
    g = A * x - b
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-10
        i = rand(1:A.n)
        x[i] = x[i] - (g[i] / A_diagonal[i])
        g = update_g(g, A, A_diagonal, i)
        total_work += (A.colptr[i+1]-A.colptr[i])
        relative_residual = norm(g) / norm(b)
    end
    return x, total_work
end

function gradient_descent(A, b)
    nonzeros_A = length(A.nzval)
    x = rand(A.n)
    g = A * x - b
    Ag = A * g
    total_work = (2 * nonzeros_A)
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-10
        alpha = dot(g, g) / (dot(g, Ag))
        x = x - alpha * g

        # Update g and Ag
        # The order of the following two statements is important
        g = g - alpha * Ag
        Ag = A * g
        relative_residual = norm(g) / norm(b)
        total_work += nonzeros_A
    end
    return x, total_work
end

# ---------------------------------- PART 2 ----------------------------------
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

# Try Gradient Descent / Coordinate Descent Algorithms
x1, total_work_1 = cyclic_coordinate_descent(T' - I, -b)
x2, total_work_2 = random_coordinate_descent(T' - I, -b)
x3, total_work_3 = gradient_descent(T' - I, -b)
println(norm(x1-x))
println(norm(x2-x))
println(norm(x3-x))
