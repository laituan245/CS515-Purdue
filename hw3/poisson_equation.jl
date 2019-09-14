using LinearAlgebra
using SparseArrays

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
uvec = A \ fvec

function run_richardson_method(A, b, alpha, iters)
    x = b
    for i=1:iters
        error = b - A * x
        x = x + alpha * error
    end
    return x
end
estimated_uvec = run_richardson_method(A, fvec, 0.249, 1000)
println(norm(estimated_uvec-uvec))


# Plotting Code
using Plots

function make_plot(uvec, n)
    N = (n+1)^2
    x = zeros(N); y = zeros(N); z = zeros(N); index = 1
    h = 1.0/(n)
    G = reshape(1:N, n+1, n+1)'
    for i=0:n
        for j=0:n
            x[index] = i * h
            y[index] = j * h
            z[index] = uvec[G[i+1,j+1]]
            index += 1
        end
    end
    return surface(x, y, z)
end

savefig(make_plot(uvec, 10), "surface.png")
