# Reference: https://www.shapeoperator.com/2014/02/22/bisecting-floats/
function bisect_step(fn::Function, x1::Float64, x2::Float64)
  xm = (x1 + x2)/2

  # Return the sub-interval with
  # oppositely-signed endpoints
  if sign(fn(x1)) != sign(fn(xm))
    return x1, xm
  else
    return xm, x2
  end
end

function bisect_root(fn::Function, x1::Float64, x2::Float64)
  @assert sign(fn(x1)) != sign(fn(x2))
  # Stop when the mean of the endpoints
  # is equal to one of the endpoints
  while x1 < (x1 + x2)/2 < x2
    x1, x2 = bisect_step(fn, x1, x2)
  end
  return x1, x2
end

function roots(a::Float64, b::Float64, c::Float64, start_range::Float64, end_range::Float64)
    function quadratic_function(x)
        return a * x * x + b * x + c
    end
    return bisect_root(quadratic_function, start_range, end_range)
end

# We can graph the quadractice function of interest to quickly identify the
# initial ranges

# In this test a = 0.2 and b = -2.6 and c = 7.7
println(roots(0.2, -2.6, 7.7, 0.0, 6.0))  # Root1 ~ 4.56351
println(roots(0.2, -2.6, 7.7, 6.0, 10.0)) # Root2 ~ 8.43649
