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
        return x, 0
    end

    rho_1 = 0; p = 0; iter = 0; mat_vec_count = 0;
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

        q = A*p; mat_vec_count += 1;
        alpha = rho / (p'*q)
        x = x + alpha * p              # update approximation vector

        r = r - alpha*q                # compute residual
        if norm(r) / bnrm2  <= tol     # check convergence
            break
        end
        rho_1 = rho
    end

    return x, mat_vec_count
end

# Neumann series-based solver
function neumann(A, b, tol)
    mat_vec_count = 0
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b) # make a copy of the right hand side
    r = b - A*x; mat_vec_count += 1; # compute the residual
    while norm(r) / bnrm2 > tol
        x .+= r
        r = b - A*x; mat_vec_count += 1;
    end
    return x, mat_vec_count
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
A = I-T'
b = b
x_true = A \ b

# Define the augmented system since A is non-symmetric
augmented_A = vcat(hcat(I, A), hcat(A', zeros(size(A))))
augmented_b = vec(vcat(b, zeros(size(b))))

# # Solve using the Conjugate gradient method
println("Conjugate gradient method")
x_cg, mat_vec_count_cg = cg(augmented_A, augmented_b, 1e-8)
x_cg = x_cg[141:280]
println(norm(A * x_cg - b) / norm(b))
println(mat_vec_count_cg)

# Solve using the neumann series-based solver
println("\nNeumann series-based solver")
x_neumann, mat_vec_count_neumann = neumann(A, b, 1e-8)
println(norm(A * x_neumann - b) / norm(b))
println(mat_vec_count_neumann)
