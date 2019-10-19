using LinearAlgebra

function solve(A, B)
    m,n = size(A)

    # Base case
    if n == 1
        return B[1] / A[1,1]
    else
        # Partition the matrix A
        k = Int32(n/2)
        A1 = A[1:k, 1:k]
        A2 = A[1:k, k+1:n]
        A3 = A[k+1:n, 1:k]
        A4 = A[k+1:n, k+1:n]

        # Partition B
        b1 = B[1:k]
        b2 = B[k+1:n]

        # Find inverse of A1
        A1_inverse = inv(A1)

        # Pre-compute C and D
        C = A1_inverse * b1
        D = A1_inverse * A2

        # Recursively solve for x2. After that, solve for x1
        x2 = solve(A4 - A3 * D, b2 - A3 * C)
        x1 = C - D * x2

        return vcat(x1,x2)
    end

end

A = rand(2, 2)
b = rand(2)
true_x = A \ b
computed_x = solve(A, b)
println(norm(true_x - computed_x))
