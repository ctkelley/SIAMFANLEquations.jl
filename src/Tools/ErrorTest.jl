#
# The functions in this file look at the status of the iteration at the
# end and set idid and errcode. 
#
# Nothing exciting in here, but it must be done.
#
"""
NewtonOK: Figure out idid and errcode for Newton's method
"""
function NewtonOK(resnorm, iline, tol, toosoon, itc, ItRules)
    maxit = ItRules.maxit
    armmax = ItRules.armmax
    printerr = ItRules.printerr
    resfail = (resnorm > tol)
    idid = ~(resfail || toosoon)
    errcode = 0
    if ~idid
        errcode =
            NewtonError(resfail, iline, resnorm, toosoon, tol, itc, maxit, armmax, printerr)
    end
    return (idid, errcode)
end

"""
PTCOK: Figure out idid and errcode
"""
function PTCOK(resnorm, tol, toosoon, ItRules, printerr)
    maxit = ItRules.maxit
    pdt0 = ItRules.pdt0
    errcode = 0
    resfail = (resnorm > tol)
    idid = ~(resfail || toosoon)
    errcode = 0
    if ~idid
        (errcode = PTCError(resnorm, maxit, pdt0, toosoon, tol, printerr))
    end
    return (idid, errcode)
end

