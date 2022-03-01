function PTCKrylovinit(
    x0,
    dx,
    F!,
    Jvec,
    delta0,
    Pvec,
    PvecKnowsdelta,
    pside,
    lsolver,
    eta,
    fixedeta,
    lmaxit,
    maxit,
    printerr,
    pdata,
)
    #
    #   Initialize the PTC-Krylov iteration.
    #
    n = length(x0)
    x = copy(x0)
    Krylov_Data = nkl_init(n, lsolver)
    kl_store = Krylov_Data.kl_store
    knl_store = Krylov_Data.knl_store
    ItRules = (
        dx = dx,
        f = F!,
        Jvec = Jvec,
        delta0 = delta0,
        Pvec = Pvec,
        PvecKnowsdelta = PvecKnowsdelta,
        pside = pside,
        lsolver = lsolver,
        kl_store = kl_store,
        knl_store = knl_store,
        eta = eta,
        fixedeta = fixedeta,
        lmaxit = lmaxit,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
    )
    return (ItRules, x, n)
end

"""
PTCUpdatei(FPS::AbstractArray, FS, x, ItRules, step, residm, delta)

Updates the PTC-Krylov iteration. This is a much simpler algorithm 
than Newton-Krylov. In particular, there is no line search to manage.

Do not mess with this function!


"""
#function PTCUpdatei(FPS::AbstractArray, FS, x, ItRules, step, residm, delta, etag)
function PTCUpdatei(FPS, FS, x, ItRules, step, residm, delta, etag)
    T = eltype(FPS)
    F! = ItRules.f
    pdata = ItRules.pdata
    #
    #
    #    step .= -(FPF \ FS)
    step .*= 0.0
    #
    #   If the preconditioner can use delta, tell it what delta is.
    #
    PvecKnowsdelta = ItRules.PvecKnowsdelta
    delta2pdata(PvecKnowsdelta, delta, pdata)
    #
    kout = Krylov_Step!(step, x, FS, FPS, ItRules, etag, delta)
    Lstats = kout.Lstats
    step = kout.step
    #
    # update solution/function value
    #
    x .= x + step
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    #
    # Update delta
    #
    delta *= (residm / resnorm)
    return (x, delta, FS, resnorm, Lstats)
end

"""
delta2pdata(PvecKnowsdelta, delta, pdata)

If your preconditioner is aware of the pseuto time step (delta)
put the value where it's supposed to be inside the precomputed data.

I also check that this has been done right and complain if not.
"""
function delta2pdata(PvecKnowsdelta, delta, pdata)
    PvecKnowsdelta || return
    # Once you're here you've told me that the preconditioner is delta-aware.
    # I will look for the array deltaval before I write to it.
    T = typeof(pdata)
    Pnames = fieldnames(T)
    valok = false
    for ip in Pnames
        valok = valok || :deltaval == ip
    end
    valok ? (pdata.deltaval[1] = delta) :
    error("PvecKnowsdelta is set to true, but there the array 
           deltaval is not a field of pdata. Check the docstrings.")
end
