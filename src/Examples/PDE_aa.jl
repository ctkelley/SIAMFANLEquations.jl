"""
PDE_aa(n=31, m=3; Vstore=Float64[], pdata=nothing, beta=1.0, maxit=40)
Solve preconditioned convection-diffusion equation with hardwired left
preconditioner using Anderson acceleration.

If you're putting this in a loop, you should allocate Vstore to 
zeros(n*n,2*(mmax+1)). Otherwise I will make a decision for you and
allocate for Vstore with each call to this function. The story on
pdata is the same. If you are calling this several times with the
same value of n, build pdata outside the call.
"""
function PDE_aa(n = 31, m = 3; Vstore = Float64[], pdata = nothing, beta = 1.0, maxit = 40)
    #
    # Process Vstore and pdata
    #
    (pdata != nothing) || (pdata = pdeinit(n))
    (length(Vstore) > 0) || (Vstore = zeros(n * n, 3 * m + 3))
    (mv, nv) = size(Vstore)
    dimvtest = ((mv == n * n) && (nv >= 2 * m + 4))
    dimvtest || error("Vstore too small")
    #
    # Call aasol and return the results.
    #
    u0 = zeros(n * n)
    rtol = 1.e-8
    atol = 1.e-8
    aout = aasol(
        hardleftFix!,
        u0,
        m,
        Vstore;
        pdata = pdata,
        maxit = maxit,
        rtol = rtol,
        atol = atol,
        beta = beta,
    )
    return aout
end
