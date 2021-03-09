function PTCKrylovinit(
    x0,
    dx,
    F!,
    Jvec,
    pdt0,
    Pvec,
    PvecKnowspdt,
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
    ItRules = (
        dx = dx,
        f = F!,
        Jvec = Jvec,
        pdt0 = pdt0,
        Pvec = Pvec,
        PvecKnowspdt = PvecKnowspdt,
        pside = pside,
        lsolver = lsolver,
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
PTCUpdatei(FPS::AbstractArray, FS, x, ItRules, step, residm, pdt)

Updates the PTC-Kyrlov iteration. This is a much simpler algorithm 
than Newton-Krylov. In particular, there is no line search to manage.

Do not mess with this function!


"""
function PTCUpdatei(FPS::AbstractArray, FS, x, ItRules, step, residm, pdt, etag)
    T = eltype(FPS)
    F! = ItRules.f
    pdata = ItRules.pdata
    #
    #
    #    step .= -(FPF \ FS)
    step .*= 0.0
#
#   If the preconditioner can use pdt, tell it what pdt is.
#
    PvecKnowspdt = ItRules.PvecKnowspdt
    pdt2pdata(PvecKnowspdt, pdt, pdata)
#
    kout = Krylov_Step!(step, x, FS, FPS, ItRules, etag, pdt)
    Lstats = kout.Lstats
    step = kout.step
    #
    # update solution/function value
    #
    x .= x + step
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    #
    # Update pdt
    #
    pdt *= (residm / resnorm)
    return (x, pdt, FS, resnorm, Lstats)
end

"""
pdt2pdata(PvecKnowspdt, pdt, pdata)

If your preconditioner is aware of the pseuto time step (pdt)
put the value where it's supposed to be inside the precomputed data.

I also check that this has been done right and complain if not.
"""
function pdt2pdata(PvecKnowspdt, pdt, pdata)
PvecKnowspdt || return
# Once you're here you've told me that the preconditioner is pdt-aware.
# I will look for the array pdtval before I write to it.
T=typeof(pdata)
Pnames = fieldnames(T)
valok = false
for ip in Pnames
   valok = valok || :pdtval == ip
end
valok ?  (pdata.pdtval[1] = pdt) :  error("PvecKnowspdt is set to true, but there the array pdtval is not a field of pdata. Check the docstrings.")
end
