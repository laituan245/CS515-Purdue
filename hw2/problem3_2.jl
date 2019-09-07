""" 4. Lookup element
Returns rho = A[i,j] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i, and j are the column indices. """
function csc_lookup(colptr, rowval, nzval, m, n, i, j)
    rho = 0
    for nzj=colptr[j]:colptr[j+1]-1 # for each entry in the j-th column
        index = rowval[nzj]
        value = nzval[nzj]
        if index == i
          rho = value
        end
    end
    return rho
end

""" 5. Lookup row
Returns x = A[i,:] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i is the row index . """
function csc_lookup_row(colptr, rowval, nzval, m, n, i)
    x = zeros(n)
    for j=1:length(colptr)-1 # for each column ...
        for nzi=colptr[j]:colptr[j+1]-1 # for each entry in the column
            row_index = rowval[nzi]
            value = nzval[nzi]
            if row_index == i
                x[j] = value
            end
        end
    end
    return x
end


# Test cases
using SparseArrays
A = sprandn(4,5,0.5) # 4x5 matrix
B = sprandn(10, 3, 0.8) # 10x3 matrix
C = sprandn(20, 20, 0.9) # 20x20 matrix

# Test the function csc_lookup
println("Test the function csc_lookup")
passed_test_1 = B[5,3] == csc_lookup(B.colptr, B.rowval, B.nzval, B.m, B.n, 5, 3)
passed_test_2 = C[8,10] == csc_lookup(C.colptr, C.rowval, C.nzval, C.m, C.n, 8, 10)
print(passed_test_1)
print(" ")
print(passed_test_2)
print("\n")
println(C[15, 12])
println(csc_lookup(C.colptr, C.rowval, C.nzval, C.m, C.n, 15, 12))

# Test the function csc_lookup_row
print("Test the function csc_lookup_row")
for i=1:4
    passed_test = A[i,:] == csc_lookup_row(A.colptr, A.rowval, A.nzval, A.m, A.n, i)
    println(passed_test)
    passed_test = C[i,:] == csc_lookup_row(C.colptr, C.rowval, C.nzval, C.m, C.n, i)
    println(passed_test)
end
