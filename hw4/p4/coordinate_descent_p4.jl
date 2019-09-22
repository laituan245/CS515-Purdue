using LinearAlgebra
using SparseArrays

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

## Poisson Equation
function laplacian(n::Integer, f::Function)
    N = (n+1)^2
    nz = (n+1) * (5*n+1)
    I = zeros(Int,nz)
    J = zeros(Int,nz)
    V = zeros(nz)
    fvec = zeros(N)
    # the transpose mirrors the row indexing we had before.
    G = reshape(1:N, n+1, n+1)' # index map, like we saw before;
    h = 1.0/(n)
    index = 1
    for i=0:n
        for j=0:n
            row = G[i+1,j+1]
            if i==0 || j == 0 || i == n || j == n
                # we are on a boudnary
                fvec[row] = 0.0
                # fill in entries in I,J,V and update index
                I[index] = row; J[index] = row; V[index]=-4; index+=1
                if j > 0 I[index] = row; J[index] = row-1; V[index]=1; index+=1 end
                if j < n I[index] = row; J[index] = row+1; V[index]=1; index+=1 end
                if i > 0 I[index] = row; J[index] = row-n-1; V[index]=1; index+=1 end
                if i < n I[index] = row; J[index] = row+n+1; V[index]=1; index+=1 end
            else
                fvec[row] = f(i*h, j*h)*h^2
                # fill in entries in I,J,V and update index
                I[index] = row; J[index] = row; V[index]=-4; index+=1
                I[index] = row; J[index] = row-1; V[index]=1; index+=1
                I[index] = row; J[index] = row+1; V[index]=1; index+=1
                I[index] = row; J[index] = row-n-1; V[index]=1; index+=1
                I[index] = row; J[index] = row+n+1; V[index]=1; index+=1
            end
        end
    end
    A = sparse(I,J,V,N,N)
    return A, fvec
end

function f(x, y)
    return 1
end

A, fvec = laplacian(10, f)
A, fvec = -A, -fvec

function cyclic_coordinate_descent(A, b)
    total_work = 0
    A_diagonal = extract_diagonal(A)
    x = rand(A.n)
    g = A * x - b
    i = 1
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-4
        x[i] = x[i] - (g[i] / A_diagonal[i])
        g = update_g(g, A, A_diagonal, i)
        total_work += (A.colptr[i+1]-A.colptr[i])
        relative_residual = norm(g) / norm(b)
        if (i == A.n) i = 1 else i += 1 end

        print(relative_residual)
        print(',')
        print(total_work)
        println()
    end
    return x
end

function random_coordinate_descent(A, b)
    total_work = 0
    A_diagonal = extract_diagonal(A)
    x = rand(A.n)
    g = A * x - b
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-4
        i = rand(1:A.n)
        x[i] = x[i] - (g[i] / A_diagonal[i])
        g = update_g(g, A, A_diagonal, i)
        total_work += (A.colptr[i+1]-A.colptr[i])
        relative_residual = norm(g) / norm(b)

        print(relative_residual)
        print(',')
        print(total_work)
        println()
    end
    return x
end

sol = random_coordinate_descent(A, fvec)

uvec = A \ fvec
println(norm(sol-uvec))
