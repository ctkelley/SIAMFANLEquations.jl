"""
NewtonError(resfail, iline, resnorm, itc, maxit, armmax)
"""
function NewtonError(resfail, iline, resnorm, itc, maxit, armmax)
if itc >= maxit
   println("Maximum iterations (maxit) of ", maxit, " exceeded")
   println("Convergence failure: residual norm too large  ",resnorm)
   println("Try increasing maxit and checking your function and 
            Jacobian for bugs.")
end
if iline 
    println("The line search failed at iteration"," ",itc)
    println("Termination with failure")
println("Current values: maxit  =  ", maxit, ", armmax = ", armmax)
end
println("Give the history array a look to see what's happening.")
println("  ")
end

"""
PTCError(resnorm, maxit, dt0)
"""
function PTCError(resnorm, maxit, dt0)
println("PTC failure; increase maxit and/or dt0")
println("Residual norm =", "  ", resnorm)
println("Current values: maxit  =  ", maxit, ",  dt0 = ", dt0)
println("Give the history array a look to see what's happening.")
println("  ")
end
