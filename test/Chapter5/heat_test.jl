function heat_test(p = 2)
    nx = (10^p) + 1
    dout = 10^(p - 1)
    na = 40
    thetal = 1.0
    thetar = 0.0
    omega = 0.9
    tau = 1.0
    Nc = 0.05
    hn_data = heat_init(nx, na, thetal, thetar, omega, tau, Nc)
    theta0 = hn_data.bcfix
    mmax = 10
    Vstore = zeros(nx, 3 * mmax + 3)
    tol = 1.e-10
    #
    # Anderson acceleration test
    #
    aout = aasol(heat_fixed!, theta0, 0, Vstore; rtol = tol, atol = tol, pdata = hn_data)
    thetabase = aout.solution
    #println(aout.history)
    test_out = thetabase[1:dout:nx]
    bench_heat = ces_heat()
    del_heat = norm(test_out - bench_heat, Inf)
    heatokaa = (del_heat < 1.e-4)
    #println(del_heat)
    for m = 1:5
        aout =
            aasol(heat_fixed!, theta0, m, Vstore; rtol = tol, atol = tol, pdata = hn_data)
        delsol = norm(aout.solution - thetabase, Inf)
        lhist = length(aout.history)
        heatmok = (delsol < 1.e-6) && (lhist == 7)
        heatokaa = heatokaa && heatmok
        #println("m=$m. solution difference = $delsol. Iterations = $lhist")
    end
    chist = aout.stats.condhist
    ahist = aout.stats.alphanorm
    heatokaa || println("aa test for heat fails")
    #
    # Newton-GMRES
    #
    return heatokaa
end

function ces_heat()
    bench_heat = [
        1.00000e+00,
        9.18027e-01,
        8.36956e-01,
        7.53557e-01,
        6.65558e-01,
        5.71475e-01,
        4.70505e-01,
        3.62437e-01,
        2.47544e-01,
        1.26449e-01,
        0.00000e+00,
    ]
    return bench_heat
end
