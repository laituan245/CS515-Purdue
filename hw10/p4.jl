# Conjugate gradient
function cg(A, b, max_it, tol)
    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    x = copy(b)
    x[:] = 0
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
        res[iter] = norm(r) / bnrm2          # check convergence
        if res[iter] <= tol
            break
        end
        rho_1 = rho

    end

    res = res[1:iter]

    return x,res
end
