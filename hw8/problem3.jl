function inner_product(x::Vector{Float64}, y::Vector{Float64})
    s = zero(Float64)
    for i=1:length(x)
        s += (x[i] * y[i])
    end
    return s
end

println(inner_product([3.1, 1.5, 0.2], [-1.0, -2.0, 3.0]))
