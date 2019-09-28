using LinearAlgebra

A = [1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1]

for i in 1:4
    for j in 1:4
        for z in 1:4
            for k in 1:4
                if i != j && i != z && i !=k && j != z && j != k && z != k
                    B = reshape([A[i, :]; A[j, :]; A[k, :]; A[z, :]], 4, 4)' # PA
                    C = I - B
                    println(eigvals(Array(C)))
                    println('\n')
                end
            end
        end
    end
end
