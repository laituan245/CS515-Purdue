function mysum(x::Vector{Float64})
  s = zero(Float64)
  for i=1:length(x)
    s += x[i]
  end
  return s
end

println(mysum([5.3, 1.25, 3.0]))
