"""
NsoliPDE(n; fixedeta=true, eta=.1)

Solve the Elliptic PDE using nsol.jl on an n x n grid. You give me
n and (optionally) sham and resdec and I return the output of nsol.
"""
function NsoliPDE(n; eta=.1, fixedeta=true, 
               rtol=1.e-7, atol=1.e-10, Pvec=Pvec2d, pside="right")
# Get some room for the residual
    u0 = zeros(n * n)
    FV = copy(u0)
    # Get the precomputed data from pdeinit
    pdata = pdeinit(n)
    # Storage for the Krylov basis
    JV = zeros(n * n, 100)
pout = nsoli(
        pdeF!,
        u0,
        FV,
        JV,
        Jvec2d;
        rtol = rtol,
        atol = atol,
        Pvec = Pvec,
        pdata = pdata,
        eta = eta,
        fixedeta = fixedeta,
        maxit = 20,
        pside = pside,
    )
return pout
end
