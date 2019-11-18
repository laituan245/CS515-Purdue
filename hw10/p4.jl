using LinearAlgebra
using DelimitedFiles
using SparseArrays

# Conjugate gradient
function cg(A, b, tol)
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b); x[:] .= 0
    r = copy(b)
    if norm(r) / bnrm2 < tol
        return x
    end

    rho_1 = 0; p = 0; iter = 0
    while true
        iter += 1
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

# Neumann series-based solver
function neumann(A, b, tol)
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b) # make a copy of the right hand side
    r = b - A*x # compute the residual
    while norm(r) / bnrm2 > tol
        x .+= r
        r = b - A*x
    end
    return x
end

# Gradient Descent
function gradient_descent(A, b, tol)
    nonzeros_A = length(A.nzval)
    x = rand(A.n)
    g = A * x - b
    Ag = A * g
    relative_residual = norm(g) / norm(b)
    while relative_residual > tol
        alpha = dot(g, g) / (dot(g, Ag))
        x = x - alpha * g

        # Update g and Ag
        # The order of the following two statements is important
        g = g - alpha * Ag
        Ag = A * g
        relative_residual = norm(g) / norm(b)
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
x_true = A \ b

# Define the augmented system since A is non-symmetric
augmented_A = vcat(hcat(I, A), hcat(A', zeros(size(A))))
augmented_b = vec(vcat(b, zeros(size(b))))

# # Solve using the Conjugate gradient method
x_cg = cg(augmented_A, augmented_b, 1e-8)[141:280]
println(norm(A * x_cg - b) / norm(b))

# # Solve using gradient descent
x_gd = gradient_descent(A, b, 1e-8)
println(norm(A * x_gd - b) / norm(b))
