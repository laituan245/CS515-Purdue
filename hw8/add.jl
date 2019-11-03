function mysum(x::Vector{Float64})
  s = zero(Float64)
  for i=1:length(x)
    s += x[i]
  end
  return s
end

function sorted_mysum(x::Vector{Float64})
  x = sort(x)
  s = zero(Float64)
  for i=1:length(x)
    s += x[i]
  end
  return s
end

println(sorted_mysum([5.3, 1.25, 3.0]))
