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

function power_method(A, b, niter)
  n = size(A,1)
  x = b
  lam = 0.0
  for i=1:niter
    y = A*x
    lam = y'*x
    x = y./norm(y)
  end
  return lam
end

# 50 x 50 instance
A = randn(50,50)
A = A' * A
b = randn(50)
true_max_eig = eigmax(A)

for niter=1:50

    S = lanczos_tri(A,b,niter)
    S = convert(Matrix, S)
    krylov_method_eig = eigmax(S)
    power_method_eig = power_method(A, b, niter)

    krylov_error = abs(krylov_method_eig-true_max_eig)
    power_method_eig = abs(power_method_eig-true_max_eig)

    println(krylov_error)
    println(power_method_eig)
end
