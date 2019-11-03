using Printf
printred(x)=print("\x1b[1m\x1b[31m"*x*"\x1b[0m")
printgreen(x)=print("\x1b[7m\x1b[32m"*x*"\x1b[0m")
printblue(x)=print("\x1b[4m\x1b[34m"*x*"\x1b[0m")

# Helper functions
function decode_Float64(x::Float64; name="")
    s = bitstring(x)
    sgn = parse(Int,s[1],base=2)
    e = parse(Int,s[2:12],base=2)
    f = parse(Int,s[13:64],base=2)
    ffrac = f/2^52

    ebias = 1023

    if length(name) > 0
        @printf("The computer representation of %s is:\n", name)
        print("  ")
        printred(s[1:1])
        printgreen(s[2:12])
        printblue(s[13:end])
        println()
        println("  which decodes to ")
    end

    println("  | sign |   exponent  |    mantissa    ")
    print("     ");printred(s[1:1]);print("     ")
    printgreen(s[2:12]);print("   ")
    printblue(s[13:end])
    println()

    fstr = (@sprintf( "%.17f", ffrac) )[3:end]

    if e == 0
        println( @sprintf("= (-1)^(%i) x 2^(%+4i)   x 0.%s = %s ", sgn, 1-ebias, fstr, string(x)) )
    elseif e == 2^11 - 1
        println( @sprintf("= (-1)^(%i) x 2^( Inf)   x 1.%s = %s ", sgn, fstr, string(x)) )
    else
        println( @sprintf("= (-1)^(%i) x 2^(%+5i)   x 1.%s = %s ", sgn, e-ebias, fstr, string(x)) )
    end
end

macro show_float(var)
    varname = string(var)
    :(decode_Float64($var;name=$varname))
end

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

# Define the problem
a = zeros(5000)
for i=1:1:5000
    a[i] = 5000 - i + 1
end
a = a .+ 0.00001

println(mysum(a))
println(sorted_mysum(a))

@show_float(mysum(a))
@show_float(sorted_mysum(a))
@show_float(12502500.05)
