using Plots
using DelimitedFiles

# Load the Yale face database of 165 pictures
A = readdlm("Yale_64.csv",',')
labels = readdlm("Yale_64_ids.csv")

# Extract even numbered images
nb_images = size(A)[1]
even_images = A[2:2:nb_images, :]
