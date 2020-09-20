"""
PTCUpdate(FPS::AbstractArray, FS, x, ItRules, step, residm, dt)

Updates the PTC iteration. This is a much simpler algorithm that Newton.
We update the Jacobian every iteration and there is no line search
to manage. 

Do not mess with this function!
In particular do not touch the line that calls PrepareJac!.
        FPF = PrepareJac!(FPS, FS, x, ItRules,dt)
PrePareJac! builds F'(u) + (1/dt)I, factors it, and sends ptcsol that
factorization.

FPF is not the same as FPS (the storage you allocate for the Jacobian)
for a reason. FPF and FPS do not have the same type, even though they
share storage. So, FPS=PrepareJac!(FPS, FS, ...) will break things.

"""
function PTCUpdate(FPS::AbstractArray, FS, x, ItRules, step, residm, dt)
    F! = ItRules.f
    pdata = ItRules.pdata
    #
    FPF = PrepareJac!(FPS, FS, x, ItRules, dt)
    #
    step .= -(FPF \ FS)
    #
    # update solution/function value
    #
    x .= x + step
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    #
    # Update dt
    #
    dt *= (residm / resnorm)
    return (x, dt, FS, resnorm)
end

"""
PTCUpdate(df::Real, fval, x, ItRules, step, residm, dt)

PTC for scalar equations. 
"""

function PTCUpdate(df::Real, fval, x, ItRules, step, resnorm, dt)
    f = ItRules.f
    dx = ItRules.dx
    fp = ItRules.fp
    df = fpeval_newton(x, f, fval, fp, dx)
    idt = 1.0 / dt
    step = -fval / (idt + df)
    x = x + step
    fval = f(x)
    # SER
    residm = resnorm
    resnorm = abs(fval)
    dt *= (residm / resnorm)
    return (x, dt, fval, resnorm)
end

#
# These functions work fine with both scalar and vector equations.
#

function PTCinit(x0, dx, F!, J!, pdata, jfact)
    #
    #   Initialize the iteration.
    #
    n = length(x0)
    x = copy(x0)
    ItRules = (dx = dx, f = F!, fp = J!, pdata = pdata, fact = jfact)
    return (ItRules, x, n)
end

function solhistinit(n, maxit, x)
#
# If you are keeping a solution history, make some room for it.
#
    solhist = zeros(n, maxit + 1)
    @views solhist[:, 1] .= x
    return solhist
end

function PTCClose(x, FS, ithist, idid, errcode, keepsolhist, solhist=[])
if keepsolhist
        sizehist = length(ithist)
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            idid = idid,
            errcode = errcode,
        )
    end
end

