using Plots
using LinearAlgebra
using DelimitedFiles

function main()
    # Load the Yale face database of 165 pictures
    A = readdlm("Yale_64.csv",',')
    labels = readdlm("Yale_64_ids.csv")
    nb_images, nb_features = size(A)

    # Nomalize each vector to unit
    for i = 1:nb_images
        A[i,:] = A[i,:] ./ max(1e-12,norm(A[i,:]));
    end

    # Let's represent each odd-numbered image in terms of the even-numbered images
    # in a least squares
    even_images = A[2:2:nb_images, :]
    X = even_images'
    for i = 1:2:nb_images
        b = A[i, :] # the current odd-numbered image
        xhat = X \ b
        reconstruction_loss = norm(b - (X * xhat))
        println(reconstruction_loss)
    end
end

main()
