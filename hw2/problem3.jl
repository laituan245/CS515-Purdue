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
        y += value*x[index]
    end
    return y
end


""" 3. Column-column inner-product
Returns rho = A[:,i]'*A[:,j] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i, and j are the column indices. """
function csc_col_col_prod(colptr, rowval, nzval, m, n, i, j)
    # Extract the i-th column
    column_i = zeros(m)
    for nzi=colptr[i]:colptr[i+1]-1 # for each entry in the i-th column
        index = rowval[nzi]
        value = nzval[nzi]
        column_i[index] = value
    end

    # Calculate rho
    rho = 0
    for nzj=colptr[j]:colptr[j+1]-1 # for each entry in the j-th column
        index = rowval[nzj]
        value = nzval[nzj]
        rho += value * column_i[index]
    end
    return rho
end

# Test cases
using SparseArrays
A = sprandn(4,5,0.5) # 4x5 matrix
B = sprandn(10, 3, 0.8) # 10x3 matrix

# Test the function csc_transpose_matvec
println("Test the function csc_transpose_matvec")
x = randn(4)
passed_test = A'*x == csc_transpose_matvec(A.colptr, A.rowval, A.nzval, A.m, A.n, x)
print("Passed test 1? ")
println(passed_test)
println(A'*x)
println(csc_transpose_matvec(A.colptr, A.rowval, A.nzval, A.m, A.n, x))

x2 = rand(10)
passed_test = B'*x2 == csc_transpose_matvec(B.colptr, B.rowval, B.nzval, B.m, B.n, x2)
print("\nPassed test 2? ")
println(passed_test)
println(B'*x2)
println(csc_transpose_matvec(B.colptr, B.rowval, B.nzval, B.m, B.n, x2))

# Test the function csc_column_projection
println("\nTest the function csc_column_projection")
passed_test = A[:, 5]'*x == csc_column_projection(A.colptr, A.rowval, A.nzval, A.m, A.n, 5, x)
print("Passed test 3? ")
println(passed_test)
println(A[:, 5]'*x)
println(csc_column_projection(A.colptr, A.rowval, A.nzval, A.m, A.n, 5, x))

passed_test = B[:,1]'*x2 == csc_column_projection(B.colptr, B.rowval, B.nzval, B.m, B.n, 1, x2)
println("\nPassed test 4? ")
println(passed_test)

# Test the function csc_col_col_prod
println("\nTest the function csc_col_col_prod")
passed_test = A[:,3]' * A[:, 4] == csc_col_col_prod(A.colptr, A.rowval, A.nzval, A.m, A.n, 3, 4)
print("Passed test 5?")
println(passed_test)

B = sprandn(10, 10, 0.8)
println(B[:,2]' * B[:,1])
println(B[:,3]' * B[:,4])
println(B[:,5]' * B[:,1])
println(csc_col_col_prod(B.colptr, B.rowval, B.nzval, B.m, B.n, 2, 1))
println(csc_col_col_prod(B.colptr, B.rowval, B.nzval, B.m, B.n, 3, 4))
println(csc_col_col_prod(B.colptr, B.rowval, B.nzval, B.m, B.n, 5, 1))


## More tests
println("\nMuch more tests")
t_matrix = sprandn(30, 35, 0.75)

x = rand(30)
passed_test_csc_transpose_matvec = (t_matrix' * x == csc_transpose_matvec(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, x))
print("passed_test_csc_transpose_matvec? ")
println(passed_test_csc_transpose_matvec)
x = rand(30)
#println(x)
passed_test_csc_transpose_matvec = (t_matrix' * x == csc_transpose_matvec(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, x))
print("passed_test_csc_transpose_matvec? ")
println(passed_test_csc_transpose_matvec)
x = rand(30)
#println(x)
passed_test_csc_transpose_matvec = (t_matrix' * x == csc_transpose_matvec(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, x))
print("passed_test_csc_transpose_matvec? ")
println(passed_test_csc_transpose_matvec)

println("\n")
x = rand(30)
passed_csc_column_projection = (t_matrix[:, 5]'*x == csc_column_projection(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, 5, x))
println("passed_csc_column_projection? ")
println(passed_csc_column_projection)
passed_csc_column_projection = (t_matrix[:, 10]'*x == csc_column_projection(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, 10, x))
println("passed_csc_column_projection? ")
println(passed_csc_column_projection)
passed_csc_column_projection = (t_matrix[:, 20]'*x == csc_column_projection(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, 20, x))
println("passed_csc_column_projection? ")
println(passed_csc_column_projection)

println("\n")
println(t_matrix[:,5]' * t_matrix[:,10])
println(csc_col_col_prod(t_matrix.colptr, t_matrix.rowval, t_matrix.nzval, t_matrix.m, t_matrix.n, 5, 10))
