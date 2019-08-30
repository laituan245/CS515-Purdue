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
