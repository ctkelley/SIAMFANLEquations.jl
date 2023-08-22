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
    FPV = zeros(n)
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
    kout3 = nsoli(
        heqf!,
        u0,
        FS,
        FPV;
        pdata = hdata,
        rtol = rtol,
        atol = atol,
        lmaxit = 40,
        eta = 0.1,
        fixedeta = true,
        lsolver = "bicgstab",
    )
    ksol = kout.solution
    dsol = dout.solution
    ksol2 = kout2.solution
    ksol3 = kout3.solution
    soltest = norm(ksol - dsol, Inf) + norm(ksol - ksol2, Inf) + norm(ksol3 - dsol, Inf)
    solpass = (soltest < 1.e-7)
    solpass || println("solpass fails")
    kfpass = (sum(kout2.stats.ikfail) == 9)
    kfpass || println("kfpass fails")
    histdiff = (dout.history - kout.history[1:8]) ./ dout.history[1]
    histpass = (norm(histdiff, Inf) < 1.e-2)
    histpass || println("histpass fails")
    histdiffb = (kout.history - kout3.history) ./ kout.history[1]
    histpassb = (norm(histdiffb, Inf) < 1.e-2)
    histpassb || println("histpassb fails")
    nkhpass = solpass && kfpass && histpass && histpassb
    return nkhpass
end
