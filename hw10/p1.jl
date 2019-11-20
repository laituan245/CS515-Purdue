using Plots
using LinearAlgebra
using SparseArrays

n = 100
on = ones(Int64,n)
A = spdiagm(-1=>-2*on[1:end-1],0=>4*on,1=>-2*on[1:end-1])
b = ones(n)

function planerot(x)
  a,b = -x[2],x[1]
  d = sqrt(b^2+a^2)
  a,b = a/d,b/d
  G = [b -a;a b]
  y = G*x
  return G,y
end

function gmres_efficient(A,b,tol,maxit)
  n = size(b)
  beta = norm(b)

  Q = zeros(length(b),maxit+1)
  Q[:,1] = b/norm(b)

  g = zeros(maxit+1)
  g[1] = beta
  Js = zeros(2,maxit) # memory to store the Js
  hist = zeros(maxit)
  H = zeros(maxit+1,maxit)

  lasti = 1
  for i = 1:maxit
    lasti = i
    v = A*Q[:,i]
    H[i+1,i] = 0 # expand k
    for j = 1:i
      H[j,i]= v'*Q[:,j]
      v = v - H[j,i]*Q[:,j]
    end

    # at this point, v = H[i+1,i]*Q[:,i+1]
    H[i+1,i] = norm(v)
    #@show v
    #@show H[i+1,i]
    Q[:,i+1] = v/H[i+1,i]

    # apply J1, ... Ji-1 to z
    for j = 1:i-1
      Jj = [Js[1,j] Js[2,j]; -Js[2,j] Js[1,j]];
      H[j:j+1,i] = Jj*H[j:j+1,i];
    end
    # create Ji to zero out the final row
    Ji,d = planerot(H[i:i+1,i])
    H[i:i+1,i] = d
    Js[1,i] = Ji[1,1]
    Js[2,i] = Ji[1,2]

    # apply Ji to g
    g[i:i+1] = Ji*g[i:i+1]
    hist[i] = abs(g[i+1]);

    if (abs(g[i+1])/beta) < tol
      break
    end
  end

  i = lasti

  # produce the solution
  y = H[1:i,1:i]\g[1:i]
  x = Q[:,1:i]*y;

  if (norm(b-A*x)/beta) < tol
    flag = 0
  elseif (hist[i]/beta) < tol
    flag = -1
  else
    flag = 1
  end

  hist = hist[1:lasti]
  return x,hist,flag
end

# Standard CG implementation
function standard_cg(A, b, tol, max_it)
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b)
    x[:] .= 0
    r = copy(b)
    error = norm(r) / bnrm2
    rho_1 = 0
    p = 0
    res = zeros(max_it)

    if error < tol
        return x,res
    end

    iter = 1
    for iter = 1:max_it
        z  = r
        rho = dot(r,z)

        if iter > 1
            beta = rho / rho_1
            p = z + beta*p
        else
            p = z
        end

        q = A*p
        alpha = rho / (p'*q)
        x = x + alpha * p                    # update approximation vector

        r = r - alpha*q                      # compute residual
        res[iter] = norm(r)                  # check convergence
        if res[iter] / bnrm2 <= tol
            break
        end
        rho_1 = rho

    end

    return x,res
end

# Using MINRES
x_minres,hist_minres,flag_minres = gmres_efficient(A,b,1.0e-8,1000)
println(norm(b-A * x_minres)/norm(b))

# Using Stanford CG
x_cg, x_res_cg = standard_cg(A, b, 1.0e-8,1000)
println(norm(b-A * x_cg)/norm(b))
println(x_res_cg[1:25])

# Simple CG
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

function simple_cg(A,b,k)
  V,T,rho = lanczos(A,b,k)
  rhs = zeros(k)
  rhs[1] = norm(b)
  y = T \ rhs # solve the tridiagonal system
  x = V[:, 1:k]*y
  return x
end

res_simple_cg = zeros(25)
for k=1:25
    x = simple_cg(A,b,k)
    res_simple_cg[k] = norm(b-A*x)
end
println(res_simple_cg)
