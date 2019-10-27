using Plots
using LinearAlgebra
using DelimitedFiles

function main()
    # Load the Yale face database of 165 pictures
    A = readdlm("Yale_64.csv",',')
    labels = readdlm("Yale_64_ids.csv")
    nb_images, nb_features = size(A)

    println(string("Number of images is ", nb_images))
    println(string("Number of features is ", nb_features))
end

main()
