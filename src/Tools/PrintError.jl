#
# The functions in this file manage error messages.
# I'm trying to give you useful hints if the iteration fails. These hints
# may become move detailed/verbose/bloviated over time.
#
"""
NewtonError(resfail, iline, resnorm, toosoon, itc, maxit, armmax, printerr)
Figure out what error message to print if the iteration fails.
"""
function NewtonError(resfail, iline, resnorm, toosoon, tol, itc, maxit, armmax, printerr)
    errcode = 0
    ~toosoon || (errcode = Lottery_Winner(resnorm, tol, printerr))
    itc < maxit || (errcode = MaxitError(resnorm, maxit, printerr))
    ~iline || (errcode = LineSearchFailure(maxit, itc, armmax, printerr))
    errcode != 0 || println("Unknown Newton error. This is not supposed to happen.")
    #
    #
    if printerr && ~toosoon
        println("Give the history array a look to see what's happening.")
        println("  ")
    end
    return errcode
end

"""
PTCError(resnorm, maxit, delta0, toosoon, tol, printerr)
"""
function PTCError(resnorm, maxit, delta0, toosoon, tol, printerr)
    ~toosoon || (errcode = Lottery_Winner(resnorm, tol, printerr))
    #if toosoon
    #errcode = Lottery_Winner(resnorm, tol, printerr)
    #else
    if printerr && ~toosoon
        println("PTC failure; increase maxit and/or delta0")
        println("Residual norm =", "  ", resnorm)
        println("Current values: maxit  =  ", maxit, ",  delta0 = ", delta0)
        println("Give the history array a look to see what's happening.")
        println("  ")
    end
    toosoon || (errcode = 10)
    return errcode
end

"""
LineSearchFailure
"""

function LineSearchFailure(maxit, itc, armmax, printerr)
    if printerr
        println("The line search failed at iteration", " ", itc)
        println("Termination with failure")
        println("Current values: maxit  =  ", maxit, ", armmax = ", armmax)
    end
    errcode = 1
    return errcode
end

"""
MaxitError
"""
function MaxitError(resnorm, maxit, printerr)
    if printerr
        println("Maximum iterations (maxit) of ", maxit, " exceeded")
        println("Convergence failure: residual norm too large  ", resnorm)
        println("Try increasing maxit and checking your function and
                 Jacobian for bugs.")
    end
    errcode = 10
    return errcode
end
"""
Lottery_Winner(resnorm, tol)
"""
function Lottery_Winner(resnorm, tol, printerr)
    if printerr
        println("Congratulations, your initial iterate met the termination criteria.")
        println("Residual norm = ", resnorm, " Tolerance = ", tol)
        println("  ")
    end
    errcode = -1
    return errcode
end


function Krylov_Error(lmaxit, ke_report)
    if ke_report == false
        println(
 "Newton-Krylov: Linear solver did not meet termination criterion at least once.
    This does not mean the nonlinear solver will fail. lmaxit= ",
            lmaxit,
        )
    end
    return true
end
