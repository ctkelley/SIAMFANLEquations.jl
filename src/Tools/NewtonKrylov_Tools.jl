function Newton_Krylov_Init(
    x0,
    dx,
    F!,
    Jvec,
    Pvec,
    pside,
    lsolver,
    eta,
    fixedeta,
    armmax,
    armfix,
    maxit,
    lmaxit,
    printerr,
    pdata,
)
    #
    #   Initialize the iteration.
    #
    eta > 0 || error("eta must be positive")
    n = length(x0)
    x = copy(x0)
    ItRules = (
        dx = dx,
        f = F!,
        Jvec = Jvec,
        Pvec = Pvec,
        pside = pside,
        lsolver = lsolver,
        eta = eta,
        fixedeta = fixedeta,
        lmaxit = lmaxit,
        armmax = armmax,
        armfix = armfix,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
    )
    return (ItRules, x, n)
end

function Krylov_Step!(step, x, FS, FPS, ItRules, etag)
    #
    # Test for too much, too soon.
    #
    lsolver = ItRules.lsolver
    lmaxit = ItRules.lmaxit
    T = eltype(FPS)
    (nk, mk) = size(FPS)
    n = length(step)
    n == nk || error("Krylov vectors wrong length")
    lmaxit < mk || error("Restarts not enabled yet")
    lsolver == "gmres" || error(lsolver, " ", "not supported")
    #T == Float64 || error("Mixed precision not supported yet. Must figure
    #out RHS scaling.")
    Jvec = ItRules.Jvec
    Pvec = ItRules.Pvec
    pdata = ItRules.pdata
    dx = ItRules.dx
    f = ItRules.f
    fixedeta = ItRules.fixedeta
    s0 = zeros(size(step))
    side = ItRules.pside
    kdata = (pdata = pdata, dx = dx, xc = x, f = f, FS = FS, Jvec = Jvec, Pvec = Pvec)
    #
    # map the Jacobian-vector and preconditioner-vector products 
    # from nsoli format to what kl_gmres wants to see
    #
    Pvecg = Pvec2
    Jvecg = Jvec2
    Pvec == nothing && (Pvecg = Pvec)
    Jvec == dirder && (Jvecg = Jvec)
    #
    #RHS=FS
    #T == Float64 || (RHS=T.(FS))
    kout = kl_gmres(s0, FS, Jvecg, FPS, etag, Pvecg; pdata = kdata, side = side)
    step = -kout.sol
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
    pdata = kdata.pdata
    atv = EvalJV(JV, v, FS, xc, pdata)
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

function EvalJV(JV, v, FS, xc, q::Nothing)
    atv = JV(v, FS, xc)
    return atv
end

function EvalJV(JV, v, FS, xc, pdata)
    atv = JV(v, FS, xc, pdata)
    return atv
end

function dirder(v, kdata)
    pdata = kdata.pdata
    dx = kdata.dx
    F = kdata.f
    FS = kdata.FS
    xc = kdata.xc
    delx = copy(xc)
    delx .= xc + dx * v
    FPP = copy(xc)
    EvalF!(F, FPP, delx, pdata)
    atv = (FPP - FS) / dx
    return atv
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
        etaflim=.5*tol/resnorm
        if etaRes <= 0.1
            etasafe = min(etamax, etaA)
        else
            etasafe = min(etamax, max(etaA, etaRes))
        end
        etag=min(etamax, max(etasafe, etaflim))
    end
    return etag
end
