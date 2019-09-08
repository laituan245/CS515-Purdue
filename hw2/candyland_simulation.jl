using LinearAlgebra
using DelimitedFiles
using SparseArrays
using Random
using Distributions

# Input the transition probability matrix T
data = readdlm("candyland-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 140,140)

function simulation(T, start_cell)
    nb_steps = 0
    current_cell = start_cell
    while current_cell != 134
        if current_cell == 5
            current_cell = 59
        elseif current_cell == 35
            current_cell = 46
        elseif current_cell == 47
            current_cell = 137
        elseif current_cell == 86
            current_cell = 138
        elseif current_cell == 117
            current_cell = 139
        else
            next_cell_distribution = Categorical(T[:,current_cell])
            next_cell = rand(next_cell_distribution, 1)[1]
            current_cell = next_cell
            nb_steps += 1
        end
    end
    return nb_steps
end

sum, count = 0, 0
start_start = 140
for i=1:1000
    global sum
    global count
    sum += simulation(T, start_start)
    count += 1
end
average = sum / count
println(average)
