using Plots
using Random
using LinearAlgebra
using DelimitedFiles

# Read and preprocess the Yale face database
function read_and_preprocess()
    # Load the Yale face database of 165 pictures
    A = readdlm("Yale_64.csv",',')
    labels = readdlm("Yale_64_ids.csv")
    nb_images, nb_features = size(A)

    # Nomalize each vector to unit
    for i = 1:nb_images
        A[i,:] = A[i,:] ./ max(1e-12,norm(A[i,:]));
    end

    return A, labels
end

# Find the value of sigma that will make face search fail for a particular image
function find_sigma(A, labels, image_idx)
    nb_images, nb_features = size(A)

    # We use all images as the database D
    D = A'

    # Extract the image of interest
    target_image = A[image_idx, :]

    # Search for gamma
    gamma = 0
    while true
        nb_fails = 0
        for i = 1:100
            preturbed_image = target_image + gamma * randn(nb_features)
            xhat = D \ preturbed_image
            if argmax(abs.(xhat)) != image_idx
                nb_fails += 1
            end
        end
        if nb_fails >= 25
            return gamma
        end
        gamma += 0.0001
    end

end

A, labels = read_and_preprocess()

# How large does sigma need to be before we won't recognize image 67 as image 67 anymore?
sigma_for_67 = find_sigma(A, labels, 67)
println(string("For image 67, sigma should be at least ", sigma_for_67))

# How large does sigma need to be before we won't recognize image 7 anymore?
sigma_for_7 = find_sigma(A, labels, 7)
println(string("For image 7, sigma should be at least ", sigma_for_7))
