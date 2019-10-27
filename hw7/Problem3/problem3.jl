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

    centered_A = center_faces(A)

    # Compute the singular value decomposition
    F = svd(centered_A')

    # Plot the leading singular vector
    u_1 = F.U[:,1]
    u_1 = reshape(u_1, 64, 64)
    savefig(heatmap(u_1, yflip=true, color=:gray), "leading_singular_vector.png")

    # Plot u_2
    u_2 = F.U[:,2]
    u_2 = reshape(u_2, 64, 64)
    savefig(heatmap(u_2, yflip=true, color=:gray), "u2.png")

    # Plot u_3
    u_3 = F.U[:,3]
    u_3 = reshape(u_3, 64, 64)
    savefig(heatmap(u_3, yflip=true, color=:gray), "u3.png")

end

function center_faces(D)
    n, k = size(D)
    e = ones(n)
    return D - (e * e' * D) / n
end

main()
