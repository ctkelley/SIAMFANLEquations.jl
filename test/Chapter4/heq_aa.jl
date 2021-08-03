function heq_aa()
    c = 0.99
    tol = 1.e-8
    fcount = [11, 10, 10, 11, 12, 12]
    anormref = [4.0, 5.4, 5.4, 5.4, 5.4, 5.4]
    condref = [1.0, 2.e2, 1.9e5, 1.9e7, 5.5e9, 6.5e10]
    n = 500
    u0 = ones(n)
    hdata = heqinit(u0, c)
    FS = zeros(n)
    FPS = zeros(n, 20)
    rtol = tol
    atol = tol
    maxit = 100
    mmax = 6
    Vstore = zeros(n, 2 * (mmax + 1))
    itrecords = zeros(6, 4)
    houtn = nsoli(
        heqf!,
        u0,
        FS,
        FPS;
        pdata = hdata,
        rtol = rtol,
        atol = atol,
        lmaxit = 10,
        eta = 0.01,
    )
    #
    # Solve H-equation with Anderson(m) for several values of m
    #
    for m = 1:6
        #houta=aasol(HeqFix!, u0, m; maxit; pdata=hdata, rtol=rtol, atol=atol)
        houta = aasol(
            HeqFix!,
            u0,
            m,
            Vstore;
            maxit = maxit,
            pdata = hdata,
            rtol = rtol,
            atol = atol,
        )
        #
        # Keep the books.
        #
        itrecords[m, 1] = hresults(houta.solution, houtn.solution)
        itrecords[m, 2] = length(houta.history)
        itrecords[m, 3] = norm(houta.stats.condhist, Inf)
        itrecords[m, 4] = norm(houta.stats.alphanorm, Inf)
    end
    #
    # Grade the results. I only use the coefficent norm and the condition
    # numbers for my own research.
    #
    solok = (norm(itrecords[:, 1], Inf) < 1.e-7)
    solok || println("Solution error in Anderson H solve")
    histok = (norm(itrecords[:, 2] - fcount, Inf) < 1.e-5)
    histok || println("History error in Anderson H solve")
    condok = (hresults(itrecords[:, 3], condref) < 1.e-1)
    condok || println("Condition error in Anderson H solve")
    normok = (hresults(itrecords[:, 4], anormref) < 1.e-1)
    normok || println("Coefficient norm error in Anderson H solve")
    #return (itrecords = itrecords)
    return solok && histok && condok && normok
end

function hresults(x, y)
    vdiff = (x - y) ./ abs.(x)
    return norm(vdiff, Inf)
end
