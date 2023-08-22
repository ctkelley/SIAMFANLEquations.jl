#
# Test results and performance for the conductive-radiative heat
# transfer problems. We compare results against column 2 of tables
# 2 and 3 in 
#
# author="C. E. Siewert and J. R. Thomas",
# title="A Computational Method for Solving a Class of Coupled
# Conductive-Radiative Heat Transfer Problems",
# journal="J. Quant. Spectrosc. Radiat. Transfer",
# year=1991,
# volume=45,
# pages="273--281"
#
#
function heat_test()
    P1ok = heat_test_examples(2, 1.0, 0.0)
    P2ok = heat_test_examples(2, 1.0, 0.5)
    return P1ok && P2ok
end
#
function heat_test_examples(p = 2, thetal = 1.0, thetar = 0.0)
    nx = (10^p) + 1
    dout = 10^(p - 1)
    na = 40
    #    thetal = 1.0
    #    thetar = 0.5
    aa_it_len = [7, 7, 7, 7, 7]
    (thetar == 0.0) || (aa_it_len = [10, 8, 7, 8, 8])
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
    aout = aasol(
        heat_fixed!,
        theta0,
        0,
        Vstore;
        maxit = 40,
        rtol = tol,
        atol = tol,
        pdata = hn_data,
    )
    thetabase = aout.solution
    test_out = thetabase[1:dout:nx]
    bench_heat = ces_heat(thetar)
    del_heat = norm(test_out - bench_heat, Inf)
    heatokaa = (del_heat < 1.e-4)
    heatokaa || println("Wrong results for xferheat: error = $del_heat")
    for m = 1:5
        aout =
            aasol(heat_fixed!, theta0, m, Vstore; rtol = tol, atol = tol, pdata = hn_data)
        delsol = norm(aout.solution - thetabase, Inf)
        lhist = length(aout.history)
        heatmok = (delsol < 1.e-6) && (lhist == aa_it_len[m])
        heatmok || println("xferheat: aa fails for m=$m and thetar=$thetar")
        heatmok || println("lhist for AA($m) = $lhist")
        heatokaa = heatokaa && heatmok
        #println("m=$m. solution difference = $delsol. Iterations = $lhist")
    end
    chist = aout.stats.condhist
    ahist = aout.stats.alphanorm
    heatokaa || println("aa test for heat fails")
    #
    # Newton-GMRES
    #
    FS = copy(theta0)
    gout = nsoli(
        FCR_heat!,
        theta0,
        FS,
        Vstore;
        pdata = hn_data,
        rtol = tol,
        atol = tol,
        dx = 1.e-5,
        eta = 0.1,
        fixedeta = false,
        lsolver = "gmres",
    )
    ndiffg = norm(gout.solution - aout.solution, Inf)
    lghist = length(gout.history)
    heatnkgok = (ndiffg < 1.e-10) && (lghist == 4)
    heatnkgok || println("xferheat: gmres fails for thetar=$thetar")
    #
    # Newton-BiCGSTAB
    #
    bout = nsoli(
        FCR_heat!,
        theta0,
        FS,
        Vstore;
        pdata = hn_data,
        rtol = tol,
        atol = tol,
        dx = 1.e-5,
        eta = 0.1,
        fixedeta = false,
        lsolver = "bicgstab",
    )
    ndiffb = norm(bout.solution - aout.solution, Inf)
    lbhist = length(bout.history)
    heatnkbok = (ndiffb < 1.e-10) && (lbhist == 4)
    heatnkbok || println("xferheat: bicgstab fails for thetar=$thetar")
    #
    return heatokaa && heatnkgok && heatnkbok
end

function ces_heat(thetar)
    if thetar == 0.0
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
    else
        bench_heat = [
            1.00000e+00,
            9.54270e-01,
            9.11008e-01,
            8.68433e-01,
            8.25127e-01,
            7.79940e-01,
            7.31936e-01,
            6.80375e-01,
            6.24709e-01,
            5.64610e-01,
            5.00000e-01,
        ]
    end
    return bench_heat
end
