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
