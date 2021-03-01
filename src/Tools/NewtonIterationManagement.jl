"""
Newtoninit: set up Newton's method
"""
function Newtoninit(
    x0,
    dx,
    F!,
    J!,
    solver,
    sham,
    armmax,
    armfix,
    resdec,
    maxit,
    printerr,
    pdata,
    jfact,
)
    #
    #   Initialize the iteration.
    #
    n = length(x0)
    x = copy(x0)
    ItRules = (
        dx = dx,
        f = F!,
        fp = J!,
        solver = solver,
        sham = sham,
        armmax = armmax,
        armfix = armfix,
        resdec = resdec,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
        fact = jfact,
    )
    return (ItRules, x, n)
end

"""
Secantinit(x0, dx, f, solver, 
         armmax, armfix, maxit, printerr, pdata, jfact)
"""
function Secantinit(x0, dx, f, solver, armmax, armfix, maxit, printerr, pdata, jfact)
    n = length(x0)
    x = copy(x0)
    ItRules = (
        f = f,
        solver = solver,
        armmax = armmax,
        armfix = armfix,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
        fact = jfact,
    )
    return (ItRules, x, n)
end

"""
NewtonOK: Figure out idid and errcode
"""
function NewtonOK(resnorm, iline, tol, toosoon, itc, ItRules)
    maxit = ItRules.maxit
    armmax = ItRules.armmax
    printerr = ItRules.printerr
    errcode = 0
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
NewtonClose: package the output of the Newton solvers
"""
function NewtonClose(x, FS, ithist, stats, idid, errcode, keepsolhist, solhist = [])
    if keepsolhist
        sizehist = length(ithist)
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
        )
    end
end

"""
CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist = [])

Collect the solution, function value, iteration stats and send them back.
"""
function CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist = [])
stats=CollectStats(ItData)
ithist=ItData.history
    if keepsolhist
        sizehist = length(ithist)
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            stats = stats,
            idid = idid,
            errcode = errcode,
        )
    end
end


