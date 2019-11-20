using Plots
using LinearAlgebra
using SparseArrays

n = 100
on = ones(Int64,n)
A = spdiagm(-1=>-2*on[1:end-1],0=>4*on,1=>-2*on[1:end-1])
b = ones(n)
