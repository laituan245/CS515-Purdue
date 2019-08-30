using Pkg; Pkg.add(["Images","FileIO","QuartzImageIO"]) # on OSX

# Load these images in Julia and convert them into matrices.
using FileIO, Images
X1 = Float64.(load("image1.png"))
X2 = Float64.(load("image2.png"))
X3 = Float64.(load("image3.png"))

# Part 2
using LinearAlgebra
part_2_answer = tr(X1+X2+X3) # this computed the trace
println("The answer to the part 2 is:")
println(part_2_answer)

# Part 4
function crude_edge_detector(nin,nout)
  Nx = reshape(1:(nin*nin), nin, nin)
  Ny = reshape(1:(nout*nout), nout, nout)
  W = zeros(nout^2,nin^2)
  for i=1:nin
    for j=1:nin
      xi = Nx[i, j]
      yj = Ny[Int(ceil(i/2)), Int(ceil(j/2))]
      is_1 = (((i % 2) + (j % 2)) % 2 == 0)
      W[yj, xi] = is_1 ? 1 : -1
    end
  end
  return W
end

W1 = crude_edge_detector(32,16)
W2 = crude_edge_detector(16,8)

println("About to start Part 5")
# Part 5
function net(x)
  global W1
  global W2
  out1 = max.(0, W1 * x)
  out2 = max.(0, W2 * out1)
  final_out = sum(out2)
  return final_out
end

@show net(vec(Float64.(load("image1.png"))))
@show net(vec(Float64.(load("image2.png"))))
@show net(vec(Float64.(load("image3.png"))))

println("\n")
println("About to start Part 7")
# Part 7
using SparseArrays
function sparse_crude_edge_detector(nin,nout)
  Nx = reshape(1:(nin*nin), nin, nin)
  Ny = reshape(1:(nout*nout), nout, nout)
  nnz = 2 * 2 * nout * nout
  I = zeros(Int, nnz) # the row index
  J = zeros(Int, nnz) # the column index
  V = zeros(nnz) # the value
  index = 1
  for i=1:nin
    for j=1:nin
      I[index] = Ny[Int(ceil(i/2)), Int(ceil(j/2))]
      J[index] = Nx[i, j]
      is_1 = (((i % 2) + (j % 2)) % 2 == 0)
      V[index] = is_1 ? 1 : -1
      index += 1
    end
  end
  return sparse(I,J,V,nout^2,nin^2)
end

# Tests
sW1 = sparse_crude_edge_detector(32,16)
W1 = crude_edge_detector(32,16)
println(sW1 == W1)

sW2 = sparse_crude_edge_detector(16,8)
W2 = crude_edge_detector(16,8)
println(sW2 == W2)
