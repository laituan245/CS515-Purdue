""" Pedagogical codes to
solve the Poisson equation using a Multigrid solver
using a Neumann boundary condition. """
# David F. Gleich
# Purdue University
# 2018-11-26

function poisson_setup(nx,ny,f)
  xs = range(0,stop=1, length=nx+1)
  ys = range(0,stop=1, length=ny+1)

  ## Create the RHS
  hx = xs[2] - xs[1]
  hy = ys[2] - ys[1]
  F = zeros(nx+2,ny+2)
  for j=2:size(F,2)-1
    for i=2:size(F,1)-1
      F[i,j] = f(xs[i], ys[j])*hx*hy
    end
  end

  return F
end

function laplacian(R)
  nx = size(R,1)
  ny = size(R,2)
  nz = 5*(nx-2)*(ny-2) + 2*nx + 2*ny - 4
  I = zeros(Int, nz)
  J = zeros(Int, nz)
  V = zeros(nz)
  map = reshape(1:nx*ny, nx, ny)
  index = 1
  b = zeros(nx*ny)
  for j=1:size(R,2)
    for i=1:size(R,1)
      row = map[i,j]
      if i==1 || i==size(R,1) || j==1 || j==size(R,2) # boundary!
        b[row] = 0.0
        I[index] = row
        J[index] = row
        V[index] = 1
        index += 1
      else
        b[row] = R[i,j]
        # fill in A[row,:]
        I[index] = row
        J[index] = row
        V[index] = 4
        index += 1

        I[index] = row
        J[index] = map[i-1,j]
        V[index] = -1
        index += 1

        I[index] = row
        J[index] = map[i+1,j]
        V[index] = -1
        index += 1

        I[index] = row
        J[index] = map[i,j-1]
        V[index] = -1
        index += 1

        I[index] = row
        J[index] = map[i,j+1]
        V[index] = -1
        index += 1
      end
    end
  end
  A = sparse(I,J,V,nx*ny,nx*ny)
  return A,b
end

function solve_poisson_direct(F)
  L,b = laplacian(F)
  return reshape(L\b, size(F)...)
end

## Now, let's see an interpolated version
""" Give the interpolated index map that preserves odd indices
1 -> 1, 1
2 -> 1, 2
3 -> 2, 2
4 -> 2, 3
...
so that X[i,j] = 1/4*X2[I1]
"""
function _interp_indexes(i)
  i,ri = divrem(i+1,2)
  return i, i+ri
end
function interpolate(X2)
  nx = (size(X2,1))+size(X2,1)-1
  ny = (size(X2,2))+size(X2,2)-1
  X = zeros(nx,ny)
  for j=1:size(X,2)
    J1,J2 = _interp_indexes(j)
    for i=1:size(X,1)
      I1,I2 = _interp_indexes(i)
      X[i,j] = (X2[I1,J1] + X2[I2,J1] + X2[I1,J2] + X2[I2,J2])/4
    end
  end
  return X
end

""" Apply the Poisson Jacobi operator.
The Jacobi iteration is D^{-1} b - D^{-1} (A-D) x^{k}
"""
function apply_poisson_jacobi(X,B)
  Y = zeros(size(X)...)
  for j=2:size(B,2)-1
    for i=2:size(B,1)-1
      Y[i,j] = B[i,j]/4 +(X[i+1,j] + X[i-1,j] + X[i,j-1] + X[i,j+1])/4
    end
  end
  return Y
end


function poisson_residual(X,B)
  R = zeros(size(X)...)
  for j=2:size(B,2)-1
    for i=2:size(B,1)-1
      R[i,j] = B[i,j] + (X[i+1,j] + X[i-1,j] + X[i,j-1] + X[i,j+1]) - 4*X[i,j]
    end
  end
  return R
end

function restrict(X)
  nyhalf = div((size(X,2) + 1) , 2)
  nxhalf = div((size(X,1) + 1) , 2)
  X2 = zeros(nxhalf,nyhalf)
  for J=3:2:size(X,2)-2
    j = div(J , 2 )
    for I=3:2:size(X,1)-2
      i = div(I , 2)
      # multiply by 4 here to accomodate the change in h spacing.
      # you could also do this on the Poisson matrix, but we have
      # to model it somewhere!
      # try removing this 4* to see how slow things go!
      X2[i+1,j+1] = 4*(1X[I-1,J-1] + 2X[I,J-1] + 1X[I+1,J-1] +
                 2X[I-1,J]   + 4X[I,J]   + 2X[I+1,J] +
                 1X[I-1,J+1] + 2X[I,J+1] + 1X[I+1,J+1])/16
    end
  end
  return X2
end

function simple_multigrid(nx,ny,niter)
  F = poisson_setup(nx,ny, (x,y) -> 1)
  X = zeros(size(F)...)
  Xfull = solve_poisson_direct(F)
  for iter=1:niter
    R = poisson_residual(X,F)
    Rhalf = restrict(R)
    X .+= interpolate(solve_poisson_direct(Rhalf))
    X = apply_poisson_jacobi(X,F)
    println(iter, " ", norm(X - Xfull)/norm(Xfull))
  end
end

function apply_poisson_multigrid(R)
  if length(R) >= 256 # we do multigrid until we get to 256x256 linear systems
    X = apply_poisson_jacobi(R,R) # one smoothing pass
    R1 = poisson_residual(X,R)
    X .+= interpolate(apply_poisson_multigrid(restrict(R1)))
    X = apply_poisson_jacobi(X,R)
    return X
  else
    return solve_poisson_direct(R)
  end
end
function poisson_multigrid_error(nx,ny,niter)
  F = poisson_setup(nx,ny, (x,y) -> 1)
  X = zeros(size(F)...)
  Xfull = solve_poisson_direct(F)
  for iter=1:niter
    R = poisson_residual(X,F)
    Rhalf = restrict(R)
    X .+= interpolate(apply_poisson_multigrid(Rhalf))
    X = apply_poisson_jacobi(X,F)
    println(iter, " ", norm(X - Xfull)/norm(Xfull), " ", norm(R)/norm(F))
  end
end

function poisson_multigrid_residual(nx,ny,niter)
  F = poisson_setup(nx,ny, (x,y) -> 1)
  X = zeros(size(F)...)
  nF = norm(F)
  for iter=1:niter
    R = poisson_residual(X,F)
    println(iter, " ", norm(R)/nF)
    if norm(R)/nF < 1.0e-8
      break
    end
    Rhalf = restrict(R)
    X .+= interpolate(apply_poisson_multigrid(Rhalf))
    X = apply_poisson_jacobi(X,F)
  end
  return X
end
