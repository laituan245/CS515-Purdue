using LinearAlgebra
using SparseArrays

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

function gradient_descent(A, b)
    nonzeros_A = length(A.nzval)
    x = rand(A.n)
    g = A * x - b
    Ag = A * g
    total_work = (2 * nonzeros_A)
    relative_residual = norm(g) / norm(b)
    while relative_residual > 1e-4
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

sol, total_work = gradient_descent(A, fvec)
uvec = A \ fvec
