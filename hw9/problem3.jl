using LinearAlgebra
using SparseArrays

function lanczos_tri(A,b,k)
  n = size(A,1)
  V = zeros(n,k+1)
  T = Tridiagonal(zeros(k-1), zeros(k), zeros(k-1))
  rho = 0.0
  vold = 0.0
  v = b/norm(b)
  for j=1:k
    y = A*v

    if j>1
      T[j-1,j] = vold'*y
      y -= T[j-1,j]*vold
    end
    T[j,j] = v'*y
    y -= T[j,j]*v

    vold = v
    rho = norm(y)
    v = y/rho
    if j < k
      T[j+1,j] = rho
    end
  end
  return T
end

A = randn(10,10)
A = A + A'
b = randn(10)
S = lanczos_tri(A,b,4)
S = convert(Matrix, S)
println(eigmax(S))
println(eigmax(A))
