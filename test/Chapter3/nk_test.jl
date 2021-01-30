"""
nk_test()

CI for nsoli
Testing Eisenstat-Walker and functions witout precomputed data
"""
function nk_test()
    passsimple = nksimple()
    passsimple || println("nksimple fails")
    jvpass = jacvec2d()
    jvpass || println("jacvec2d fails")
    nkpass = passsimple && jvpass
    return nkpass
end
"""
nksimple()

Test nsoli with the simple 2D problem and line search failure and success.
"""
function nksimple()
    x0 = [2.0; 0.5]
    FPS = zeros(2, 4)
    FPJ = zeros(2, 2)
    FS = copy(x0)
    #
    # For the easy problem we will do analytic Jacobians for
    # Newton and forward difference directional derivatives for Newton-GMRES
    #
    dout = nsol(simple!, x0, FS, FPJ, jsimple!; sham = 1, keepsolhist = true)
    kout = nsoli(
        simple!,
        x0,
        FS,
        FPS;
        eta = 1.e-10,
        lmaxit = 3,
        keepsolhist = true,
        fixedeta = false,
    )
    dsolhist = norm(kout.solhist - dout.solhist, Inf)
    shpass = (dsolhist < 1.e-7)
    shpass || println("solhist compare fails in easy nksimple", dsolhist)
    vconverge = krstest(dout, kout)
    #
    # For the stagnating problem we will do analytic Jacobians for
    # Newton and analytic Jacobian-vector products for Newton-GMRES
    #
    x0 = [3.0; 5.0]
    dout = nsol(simple!, x0, FS, FPJ, jsimple!; sham = 1)
    kout = nsoli(simple!, x0, FS, FPS, Jvsimple; fixedeta = true, eta = 1.e-10, lmaxit = 2)
    vdiverge = krstest(dout, kout)
    vdiverge || println("failure hard nksimple problem")
    return vconverge && vdiverge && shpass
end

function krstest(dout, kout)
    hdiff = norm(kout.history - dout.history, Inf)
    hpass = (hdiff < 5.e-7)
    hpass || println("history compare fails in nksimple")
    #
    adiff = kout.stats.iarm - dout.stats.iarm
    apass = (sum(adiff) == 0)
    apass || println("line search compare fails in nksimple")
    #
    fdiff = kout.stats.ifun - dout.stats.ifun
    fpass = (sum(fdiff) == 0)
    fpass || println("function value compare fails in nksimple")
    #
    soldiff = kout.solution - dout.solution
    solpass = (norm(soldiff, Inf) < 1.e-9)
    solpass || println("solution compare fails in nksimple  ", norm(soldiff, Inf))
    krpass = (fpass && apass && hpass && solpass)
end

function Jvsimple(v, FS, x)
    #println(typeof(v))
    JM = zeros(2, 2)
    jsimple!(JM, FS, x)
    Jvec = JM * v
    return Jvec
end


"""
jacvec2d()

Analytic Jacobian-vector product. Compare Eisenstat-Walker to 
fixed eta. Test precomputed data support.
"""
function jacvec2d()
    x0 = ones(2)
    fv = zeros(2)
    jv = zeros(2, 2)
    jvs = zeros(2, 3)
    pdata = zeros(2)
    nout = nsol(f!, x0, fv, jv; sham = 1, pdata = pdata)
    kout = nsoli(f!, x0, fv, jvs, JVec; fixedeta = false, eta = 0.9, lmaxit = 2, pdata = pdata)
    kout2 = nsoli(fv2!, x0, fv, jvs, JVecv2; fixedeta = true, eta = 0.1, 
                 lmaxit = 2)
    histdiff = norm(nout.history - kout2.history)
    histpass = (histdiff < 1.e-5)
    histpass || println("hist test fails in jacvec2d")
    ncost = funcost(nout)
    nplot = acost(nout)
    kcost = funcost(kout)
    kplot = acost(kout)
    kcost2 = funcost(kout2)
    kplot2 = acost(kout2)
    costpass = (ncost == 10) && (kcost == 15) && (kcost2 == 14)
    costpass || println("cost compare fails in jacvec2d")
    costpass || println(ncost,"  ", kcost, "  ", kcost2)
    soldiff = (
        norm(kout.solution - nout.solution, Inf) +
        norm(kout2.solution - nout.solution, Inf)
    )
    solpass = (soldiff < 1.e-7)
    solpass || println("solution compare fails in jacvec2d")
    jvpass = histpass && costpass && solpass
    return jvpass
end

function f!(fv, x, pdata)
    fv[1] = x[1] + sin(x[2])
    fv[2] = cos(x[1] + x[2])
end


"""
fv2!(fv, x)

Function evaluation witout precomputed data for testing.
"""
function fv2!(fv, x)
    fv[1] = x[1] + sin(x[2])
    fv[2] = cos(x[1] + x[2])
end

"""
JVecv2(v, fv, x)

Precondition without precomputed/stored data
"""
function JVec(v, fv, x)
    jvec = zeros(2)
    p = -sin(x[1] + x[2])
    jvec[1] = v[1] + cos(x[2]) * v[2]
    jvec[2] = p * (v[1] + v[2])
    return pdata
end


"""
JVec(v, fv, x, pdata)

Precondition with precomputed/stored data
"""
function JVec(v, fv, x, pdata)
    jvec = zeros(2)
    p = -sin(x[1] + x[2])
    pdata[1] = v[1] + cos(x[2]) * v[2]
    pdata[2] = p * (v[1] + v[2])
    return pdata
end

function funcost(itout)
    netcost = itout.stats.ifun + itout.stats.ijac
    cost = sum(netcost)
end

function acost(itout)
    netcost = itout.stats.ifun + itout.stats.ijac
    cost = cumsum(netcost)
end
