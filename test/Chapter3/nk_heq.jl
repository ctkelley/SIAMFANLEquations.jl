"""
nk_heq()

CI for nsoli and H-equation.
"""
function nk_heq()
    n = 32
    u0 = zeros(n)
    FS = zeros(n)
    FPS = zeros(n, 20)
    FPJ = zeros(n, n)
    c = 0.999
    atol = 1.e-9
    rtol = 1.e-9
    hdata = heqinit(u0, c)
    dout =
        nsol(heqf!, u0, FS, FPJ, heqJ!; rtol = rtol, atol = atol, pdata = hdata, sham = 1)
    kout = nsoli(
        heqf!,
        u0,
        FS,
        FPS;
        pdata = hdata,
        rtol = rtol,
        atol = atol,
        lmaxit = -1,
        eta = 0.1,
        fixedeta = false,
    )
    kout2 = nsoli(
        heqf!,
        u0,
        FS,
        FPS;
        pdata = hdata,
        rtol = rtol,
        atol = atol,
        lmaxit = 2,
        eta = 0.01,
    )
    ksol = kout.solution
    dsol = dout.solution
    ksol2 = kout2.solution
    soltest = norm(ksol - dsol, Inf) + norm(ksol - ksol2, Inf)
    solpass = (soltest < 1.e-7)
    kfpass = (sum(kout2.stats.ikfail) == 9)
    histdiff = (dout.history - kout.history[1:8]) ./ dout.history[1]
    histpass = (norm(histdiff, Inf) < 1.e-2)
    nkhpass = solpass && kfpass && histpass
    return nkhpass
end
