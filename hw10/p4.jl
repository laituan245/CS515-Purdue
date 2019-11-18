using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Conjugate gradient
function cg(A, b, tol)
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b); x[:] = 0
    r = copy(b)
    if norm(r) / bnrm2 < tol
        return x
    end

    rho_1 = 0; p = 0
    while true
        z = r
        rho = dot(r,z)

        if iter > 1
            beta = rho / rho_1
            p = z + beta*p
        else
            p = z
        end

        q = A*p
        alpha = rho / (p'*q)
        x = x + alpha * p              # update approximation vector

        r = r - alpha*q                # compute residual
        if norm(r) / bnrm2  <= tol     # check convergence
            break
        end
        rho_1 = rho
    end

    return x
end

# Build the Candyland linear system
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

# Prepare the vector b
b = ones((140, 1))
b[134] = 0

# Define the system
A = T' - I
b = -b
