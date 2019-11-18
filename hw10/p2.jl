using Plots
using LinearAlgebra
using SparseArrays

# Define constants
lambda_1 = 0.1
lambda_n = 100
p = 0.9
n = 30

# Function for defining a Strako\v{s} matrix
function build_strakos_matrix(lambda_1, lambda_n, p, n)
    d = zeros(n)
    for i = 1:n
        d[i] = lambda_1 + ((i-1)/(n-1)) * (lambda_n - lambda_1) * p^(n-i)
    end
    return Diagonal(d)
end

A = build_strakos_matrix(lambda_1, lambda_n, p, n)

# Function for the Lanczos method
function lanczos(A,b,k)
  n = size(A,1)
  V = zeros(n,k+1)
  T = Tridiagonal(zeros(k-1), zeros(k), zeros(k-1))
  rho = 0.0
  V[:,1] = b/norm(b)
  for j=1:k
    y = A*V[:,j]
    for i=max(j-1,1):j
      T[i,j] = V[:,i]'*y
      y -= T[i,j]*V[:,i]
    end
    rho = norm(y)
    V[:,j+1] = y/rho
    if j < k
      T[j+1,j] = rho
    end
  end
  return V,T,rho
end

# Part 1
quantities = zeros(n)
b = ones(n) / sqrt(n)
for k=1:n
    V = lanczos(A, b, k)[1][:,1:k]
    quantities[k] = log10(norm(V' * V - I) + 1e-20)
end
savefig(plot(quantities, linewidth=2),"part1.png")

# Part 4
function arnoldi(A,b,k)
  n = size(A,1)
  V = zeros(n,k+1)
  H = zeros(k+1,k)
  V[:,1] = b/norm(b)
  for j=1:k
    y = A*V[:,j]
    for i=1:k
      H[i,j] = V[:,i]'*y
      y -= H[i,j]*V[:,i]
    end
    H[j+1,j] = norm(y)
    V[:,j+1] = y/H[j+1,j]
  end
  return V,H
end

quantities = zeros(60)
b = rand(n); b = b / norm(b)
for k=1:60
    V, T = arnoldi(A, b, k)
    quantities[k] = log10(norm(A * V[:,1:k] - V * T) + 1e-20)
end
savefig(plot(quantities, linewidth=2),"part4.png")
