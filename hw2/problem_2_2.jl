function laplacian(n::Integer, f::Function)
    N = (n+1)^2
    nz = <fill-in>
    I = zeros(Int,nz)
    J = zeros(Int,nz)
    V = zeros(nz)
    fvec = zeros(N)
    # the transpose mirrors the row indexing we had before.
    G = reshape(1:N, n+1, n+1)' # index map, like we saw before;
    h = 1.0/(n)
    index = 1
    for i=0:n
        for j=0:n
            row = G[i+1,j+1]
            if i==0 || j == 0 || i == n || j == n
                # we are on a boudnary
                fvec[row] = 0.0
                # fill in entries in I,J,V and update index
            else
                fvec[row] = f(i*h, j*h)*h^2
                # fill in entries in I,J,V and update index
            end
        end
    end
    A = sparse(I,J,V,N,N)
    return A, fvec
end
