"""
ptcsolsc_test

Make sure ptcsolsc finds the stable steady state sqrt(2)/2

f(u) = u^3 - lambda u = 0 with lambda = 1/2

The answer(s!!) are u=0, -sqrt(2)/2, and sqrt(2)/2.

The initial iterate is u0=.1, so the correct answer is sqrt(2)/2.

I'm testing for correctness and a match of the iterations statistics

"""
function ptcsolsc_test()
    u0 = 0.1
    ustable = 0.5 * sqrt(2.0)
    uunstable = 0.0
    lambda = 0.5
    #
    # Convergence to the right solution
    #
    ptcdata1 = ptcsolsc(sptest, u0, sptestp; dt0 = 1.0, rtol = 1.e-12)
    ptcdata2 = ptcsolsc( spitchfork, u0; dt0 = 1.0, rtol = 1.e-12, 
        pdata = lambda, keepsolhist = false,)
    ptcdatasec = secant( spitchfork, u0; rtol = 1.e-12, 
        pdata = lambda, keepsolhist = false,)
    dh = ptcdata1.history - ptcdata2.history
    ndh = norm(dh[:, 1], Inf)
    fdok = (ndh < 1.e-7)
    ptcerr = ptcdata1.solhist .- ustable
    ptcfun = ptcdata1.history
    secok = abs(ptcdatasec.solution) < 1.e-10
    solok = (abs(ptcdata1.solution - ustable) < 1.e-10)
    funok = (abs(ptcdata1.functionval) < 1.e-12)
    histok = (length(ptcfun) == 18)
    ptcdataf = ptcsolsc(sptest, u0; dt0 = 0.1, rtol = 1.e-12)
    errcode = ptcdataf.errcode
    failok = ~ptcdataf.idid && (errcode == 10)
    ptcok = fdok && solok && funok && histok && failok && secok
    if ~ptcok
        println("Failure in Scalar PTC")
        println(ptcdata1)
    end
    return ptcok
end
