using LinearAlgebra
using SparseArrays

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
