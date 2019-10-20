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
    D = even_images'
    min_reconstruction_loss, min_idx = Inf, nothing
    max_reconstruction_loss, max_idx = -Inf, nothing
    for i = 1:2:nb_images
        b = A[i, :] # the current odd-numbered image
        xhat = D \ b
        reconstruction_loss = norm(b - (D * xhat))
        if reconstruction_loss < min_reconstruction_loss
            min_reconstruction_loss = reconstruction_loss
            min_idx = i
        end
        if reconstruction_loss > max_reconstruction_loss
            max_reconstruction_loss = reconstruction_loss
            max_idx = i
        end
    end

    println(min_idx)
    println(max_idx)
end

main()
