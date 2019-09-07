""" 1. Sparse matrix-transpose multiplication by a vector
Returns y = A'*x where A is given by the CSC arrays
colptr, rowval, nzval, m, n and x is the vector. """
function csc_transpose_matvec(colptr, rowval, nzval, m, n, x)
  y = zeros(n)
  for j=1:length(colptr)-1 # for each column ...
      for nzi=colptr[j]:colptr[j+1]-1 # for each entry in the column
          i = rowval[nzi]
          v = nzval[nzi]
          y[j] += v*x[i]
      end
  end
  return y
end


""" 2. Column inner-product
Returns  = A[:,i]'*x where A is given by the CSC arrays
colptr, rowval, nzval, m, n and x is the vector. """
function csc_column_projection(colptr, rowval, nzval, m, n, i, x)
    y = 0
    for nzi=colptr[i]:colptr[i+1]-1 # for each entry in the i-th column
        index = rowval[nzi]
        value = nzval[nzi]
        y += v*x[index]
    end
    return y
end

# Test cases
using SparseArrays
A = sprandn(4,5,0.5) # 4x5 matrix
B = sprandn(10, 3, 0.8) # 10x3 matrix

# Test the function csc_transpose_matvec
x = randn(4)
passed_test = A'*x == csc_transpose_matvec(A.colptr, A.rowval, A.nzval, A.m, A.n, x)
print("Passed test 1? ")
println(passed_test)
println(A'*x)
println(csc_transpose_matvec(A.colptr, A.rowval, A.nzval, A.m, A.n, x))

x2 = rand(10)
passed_test = B'*x2 == csc_transpose_matvec(B.colptr, B.rowval, B.nzval, B.m, B.n, x2)
print("\nPassed test2? ")
println(passed_test)
println(B'*x2)
println(csc_transpose_matvec(B.colptr, B.rowval, B.nzval, B.m, B.n, x2))
