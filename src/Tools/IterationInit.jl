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
    keepsolhist
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
        fact = jfact
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
PTCinit(x0, dx, F!, J!, pdt0, maxit, pdata, jfact, keepsolhist)

PTCinit: get organized for PTC 
"""
function PTCinit(x0, dx, F!, J!, pdt0, maxit, pdata, jfact, keepsolhist,
                  jknowsdt=false)
    #
    #   Initialize the iteration.
    #
    n = length(x0)
    x = copy(x0)
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    ItRules =
        (dx = dx, f = F!, fp = J!, pdt0 = pdt0, maxit = maxit, 
          pdata = pdata, fact = jfact, jknowsdt=jknowsdt)
    return (ItRules, x, n, solhist)
end

"""
Newton_Krylov_Init( x0, dx, F!, Jvec, Pvec, pside, lsolver, eta,
    fixedeta, armmax, armfix, maxit, lmaxit, printerr, pdata, keepsolhist)

Newton_Krylov_Init: set up nsoli
"""
function Newton_Krylov_Init( x0, dx, F!, Jvec, Pvec, pside, lsolver, eta,
    fixedeta, armmax, armfix, maxit, lmaxit, printerr, pdata, keepsolhist)
    #
    #   Initialize the iteration.
    #
    eta > 0 || error("eta must be positive")
    n = length(x0)
    x = copy(x0)
    tmp1=zeros(n,); tmp2=zeros(n,); tmp3=zeros(n,); tmp4=zeros(n,);
    kl_store = (tmp1, tmp2, tmp3, tmp4)
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    ItRules = (
        dx = dx,
        f = F!,
        Jvec = Jvec,
        Pvec = Pvec,
        pside = pside,
        lsolver = lsolver,
        kl_store = kl_store,
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

