#
# The functions in this file manage the Newton-Krylov step and
# the Jacobian/preconditioner - vector products. 
#
"""
Krylov_Step!(step, x, FS, FPS, ItRules, etag, pdt = 0)

Take a Newton-Krylov step. This function does lots of its work
mapping nonlinear problems to linear solvers. Only then do I get
to deploy the Krylov linear solvers.
"""
function Krylov_Step!(step, x, FS, FPS, ItRules, etag, pdt = 0)
    #
    # Test for too much, too soon.
    #
    lsolver = ItRules.lsolver
    lmaxit = ItRules.lmaxit
    T = eltype(FPS)
    (nk, mk) = size(FPS)
    n = length(step)
    n == nk || error("Krylov vectors wrong length")
    lsolver == "gmres" || error(lsolver, " ", "not supported")
    Jvec = ItRules.Jvec
    Pvec = ItRules.Pvec
    kl_store = ItRules.kl_store
    pdata = ItRules.pdata
    dx = ItRules.dx
    f = ItRules.f
    fixedeta = ItRules.fixedeta
    #    s0 = zeros(size(step))
    #
    #   Initial iterate for step is zero.
    #
    s0 = step
    s0 .*= 0.0
    side = ItRules.pside
    #
    # map the Jacobian-vector and preconditioner-vector products 
    # from nsoli format to what kl_gmres wants to see
    #
    kdata = (
        pdata = pdata,
        dx = dx,
        xc = x,
        f = f,
        FS = FS,
        Jvec = Jvec,
        Pvec = Pvec,
        pdt = pdt,
    )
    Pvecg = Pvec2
    Jvecg = Jvec2
    Pvec == nothing && (Pvecg = Pvec)
    Jvec == dirder && (Jvecg = Jvec)
    #
    #RHS=FS
    #T == Float64 || (RHS=T.(FS))
    kout = kl_gmres(
        s0,
        FS,
        Jvecg,
        FPS,
        etag,
        Pvecg;
        kl_store = kl_store,
        pdata = kdata,
        side = side,
        lmaxit = lmaxit,
    )
    step .= -kout.sol
    reshist = kout.reshist
    lits = kout.lits
    idid = kout.idid
    Lstats = (reshist = reshist, lits = lits, idid = idid)
    return (step = step, Lstats = Lstats)
end

function Pvec2(v, kdata)
    F = kdata.f
    FS = kdata.FS
    xc = kdata.xc
    PV = kdata.Pvec
    pdata = kdata.pdata
    ptv = EvalPV(PV, v, xc, pdata)
    return ptv
end


function Jvec2(v, kdata)
    F = kdata.f
    FS = kdata.FS
    xc = kdata.xc
    JV = kdata.Jvec
    pdt = kdata.pdt
    pdata = kdata.pdata
    atv = EvalJV(JV, v, FS, xc, pdt, pdata)
    return atv
end

function EvalPV(PV, v, xc, q::Nothing)
    ptv = PV(v, xc)
    return ptv
end

function EvalPV(PV, v, xc, pdata)
    ptv = PV(v, xc, pdata)
    return ptv
end

function EvalJV(JV, v, FS, xc, pdt, q::Nothing)
    atv = JV(v, FS, xc)
    ptcmv!(atv, v, pdt)
    #    if pdt > 0
    #        atv .= atv + (1.0 / pdt) * v
    #    end
    return atv
end

function EvalJV(JV, v, FS, xc, pdt, pdata)
    atv = JV(v, FS, xc, pdata)
    ptcmv!(atv, v, pdt)
    #    if pdt > 0
    #        atv .= atv + (1.0 / pdt) * v
    #    end
    return atv
end

function dirder(v, kdata)
    pdata = kdata.pdata
    dx = kdata.dx
    F = kdata.f
    FS = kdata.FS
    xc = kdata.xc
    pdt = kdata.pdt
    delx = copy(xc)
    delx .= xc + dx * v
    FPP = copy(xc)
    EvalF!(F, FPP, delx, pdata)
    atv = (FPP - FS) / dx
    ptcmv!(atv, v, pdt)
    #    if pdt > 0
    #        atv .= atv + (1.0 / pdt) * v
    #    end
    return atv
end

function ptcmv!(atv, v, pdt)
    (pdt == 0.0) || (atv .= atv + (1.0 / pdt) * v)
    #return atv
end

"""
forcing(itc, residratio, etag, ItRules, tol)

Compute the Eisenstat-Walker forcing term
"""
function forcing(itc, residratio, etag, ItRules, tol, resnorm)
    gamma = 0.9
    etamax = ItRules.eta
    fixedeta = ItRules.fixedeta
    if fixedeta || (itc == 0)
        etag = etamax
    else
        etaRes = gamma * (etag^2)
        etaA = gamma * (residratio^2)
        etaflim = 0.5 * tol / resnorm
        if etaRes <= 0.1
            etasafe = min(etamax, etaA)
        else
            etasafe = min(etamax, max(etaA, etaRes))
        end
        etag = min(etamax, max(etasafe, etaflim))
    end
    return etag
end
