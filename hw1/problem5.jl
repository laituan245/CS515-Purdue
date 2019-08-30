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

# Part 5
function net(x)
  global W1
  global W2
  out1 = max.(0, W1 * vec(x))
  out2 = max.(0, W2 * out1)
  final_out = sum(out2)
  return final_out
end

println(net(X1))
