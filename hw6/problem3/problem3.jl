using Plots
using DelimitedFiles

# Load the Yale face database of 165 pictures
A = readdlm("Yale_64.csv",',')
labels = readdlm("Yale_64_ids.csv")
