"""
NewtonError(resfail, iline, resnorm, itc, maxit, armmax,printerr)
"""
function NewtonError(resfail, iline, resnorm, itc, maxit, armmax,printerr)
errcode=0
if itc >= maxit
   if printerr
   println("Maximum iterations (maxit) of ", maxit, " exceeded")
   println("Convergence failure: residual norm too large  ",resnorm)
   println("Try increasing maxit and checking your function and 
            Jacobian for bugs.")
   end 
   errcode+=10
end
if iline 
   if printerr
    println("The line search failed at iteration"," ",itc)
    println("Termination with failure")
    println("Current values: maxit  =  ", maxit, ", armmax = ", armmax)
    end
    errcode+=1
end
if printerr
println("Give the history array a look to see what's happening.")
println("  ")
end
return errcode
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

"""
Lottery_Winner(resnorm, tol)
"""
#function Lottery_Winner(resnorm, tol, printerr)
#if printerr
#println("Congratulations, your initial iterate met the teremination criteria.")
#println("Residual norm = ", resnorm, " Tolerance = ", tol)
#println("  ")
#end
#errcode = -1
#return errcode
#end
