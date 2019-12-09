# Midterm: https://www.cs.cornell.edu/courses/cs6210/2019fa/hw/mt.pdf

using LinearAlgebra

# Problem 1.3
Z = rand(7, 3)
A = I + Z * Z'
println(string("Condition number of A is ", cond(A, 2)))
U,Sig,V = svd(Z, full=true)
println(string("My predicted condition number of A is ", 1.0 + Sig[1]^2))
