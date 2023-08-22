#
# Test results and performance for the conductive-radiative heat
# transfer problems. 
# 
# This test makes sure the traps and error codes for failure
# do what I want.
#
function heat_test2()
    P1ok = heat_test2_examples()
    return P1ok
end
#
function heat_test2_examples(p = 2, thetal = 1.0, thetar = 2.0, omega = 0.5, tau = 4.0)
    nx = (10^p) + 1
    na = 40
    Nc = 0.05
    hn_data = heat_init(nx, na, thetal, thetar, omega, tau, Nc)
    theta0 = hn_data.bcfix
    mmax = 50
    Vstore = zeros(nx, 3 * mmax + 3)
    tol = 1.e-10
    errcodes = [-2, 0, 10]
    errtarget = [1.e4, 1.e-8, 1.e-5]
    Pok = true
    #
    # Newton-GMRES to obtain a converged result
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
    thetabase = gout.solution
    gmhistok = (length(gout.history) == 7)
    gmjacok = (sum(gout.stats.ijac) == 28)
    gmconvok = (gout.errcode == 0)
    gmresok = gmhistok && gmjacok && gmconvok
    gmresok || println("nsoli fails in heat_test2")
    Pok = Pok && gmresok
    #
    # Anderson acceleration test
    #
    iec = 1
    for m in [5, 10, 20]
        aout = aasol(
            heat_fixed!,
            theta0,
            m,
            Vstore;
            rtol = tol,
            atol = tol,
            pdata = hn_data,
            maxit = 50,
        )
        delsol = norm(aout.solution - thetabase, Inf)
        errc = aout.errcode
        ecodeok = (aout.errcode == errcodes[iec])
        ecodeok || println("ecode test fails, heat_test2, m=$m")
        Pok = Pok && ecodeok
        delok = (delsol < errtarget[iec])
        delok || println("sol err test fails, heat_test2, m=$m")
        Pok = Pok && delok
        iec += 1
        println("For m=$m: error=$delsol, errcode = $errc")
    end
    return Pok
end
