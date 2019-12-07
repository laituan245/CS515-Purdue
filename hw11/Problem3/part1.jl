using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Random

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

# Read the eigenfaces matrix and compute svd using Julia built-in function
A = readdlm(download("http://www.cs.purdue.edu/homes/dgleich/cs515-2019/homeworks/Yale_64.csv"),',')
M = A'

# Build bipartite matrix
m, n = size(M)
B = [spzeros(m, m) M; M' spzeros(n,n)]
for k=5:16
    global(B); global(m); global(n)
    println(string("\nk = ", k))
    V, T, rho = lanczos(B,ones(m+n),k)
    display(sort(eigvals(Matrix(T)), rev=true)[1:5])
end


# Compare to the built-in SVD function of Julia
println("\nUse the built-in SVD function of Julia")
U,Sig,V = svd(M)
println(Sig[1:5])
