#
# The functions in this file initialize the iterations 
#
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
    keepsolhist,
)
    #
    #   Initialize the iteration.
    #
    n = length(x0)
    x = copy(x0)
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
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
    return (ItRules, x, n, solhist)
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
PTCinit(x0, dx, F!, J!, delta0, maxit, pdata, jfact, keepsolhist)

PTCinit: get organized for PTC 
"""
function PTCinit(x0, dx, F!, J!, delta0, maxit, pdata, jfact, keepsolhist, jknowsdt = false)
    #
    #   Initialize the iteration.
    #
    n = length(x0)
    x = copy(x0)
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    ItRules = (
        dx = dx,
        f = F!,
        fp = J!,
        delta0 = delta0,
        maxit = maxit,
        pdata = pdata,
        fact = jfact,
        jknowsdt = jknowsdt,
    )
    return (ItRules, x, n, solhist)
end

"""
Newton_Krylov_Init( x0, dx, F!, Jvec, Pvec, pside, lsolver, eta,
    fixedeta, armmax, armfix, maxit, lmaxit, printerr, pdata, u
    Krylov_Data, keepsolhist)

Newton_Krylov_Init: set up nsoli
"""
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
    Krylov_Data,
    keepsolhist,
)
    #
    #   Initialize the iteration.
    #
    eta > 0 || error("eta must be positive")
    n = length(x0)
    #
    # Not for tourists! You have the opportunity, which you should decline,
    # to allocate the internal space for gmres in the call to nsoli. Only
    # do this for continuation or IVP integration, if at all. You can break
    # stuff with this.
    #
    if Krylov_Data == nothing
        Krylov_Data = nkl_init(n, lsolver)
        #    kl_store = kstore(n,lsolver)
        #    knl_store = knlstore(n)
    end
    kl_store = Krylov_Data.kl_store
    knl_store = Krylov_Data.knl_store
    x = knl_store.xval
    x .= x0
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    ((lmaxit == -1) && (lsolver == "bicgstab")) && (lmaxit = 5)
    ItRules = (
        dx = dx,
        f = F!,
        Jvec = Jvec,
        Pvec = Pvec,
        pside = pside,
        lsolver = lsolver,
        kl_store = kl_store,
        knl_store = knl_store,
        eta = eta,
        fixedeta = fixedeta,
        lmaxit = lmaxit,
        armmax = armmax,
        armfix = armfix,
        maxit = maxit,
        printerr = printerr,
        pdata = pdata,
    )
    return (ItRules, x, n, solhist)
end



"""
solhistinit(n, maxit, x)

Am I keeping the solution history? If so, allocate the space.
"""
function solhistinit(n, maxit, x)
    #
    # If you are keeping a solution history, make some room for it.
    #
    solhist = zeros(n, maxit + 1)
    @views solhist[:, 1] .= x
    return solhist
end
