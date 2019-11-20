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

# Using MINRES
x_minres,hist_minres,flag_minres = gmres_efficient(A,b,1.0e-8,1000)
println(norm(b-A * x_minres)/norm(b))
