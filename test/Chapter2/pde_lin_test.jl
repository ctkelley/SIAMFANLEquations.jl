"""
pde_lin_test(n)

Test the linear operators.
"""
function pde_lin_test(n)
    h = 1.0 / (n + 1)
    x = collect(h:h:1.0-h)
    fdata = fishinit(n)
    z = rand(n, n)
    L2d = Lap2d(n)
    DX = Dx2d(n)
    DY = Dy2d(n)
    lapok = lap_test(n, x, z, L2d, fdata)
    discok = disc_test(n, x, L2d, DX, DY, fdata)
    lapok && discok
end

"""
disc_test()
Have I broken the discretizations?
"""
function disc_test(n, x, L2d, DX, DY, fdata)
    n2 = n * n
    ue = solexact(x)
    ux = dxexact(x)
    uy = dyexact(x)
    D2u = l2dexact(x)
    ue1 = reshape(ue, (n2,))
    uex1 = reshape(ux, (n2,))
    uey1 = reshape(uy, (n2,))
    ued21 = reshape(D2u, (n2,))
    # test DX
    dx21 = DX * ue1
    dxerr = norm(dx21 - uex1, Inf)
    # test DY
    dy21 = DY * ue1
    dyerr = norm(dy21 - uey1, Inf)
    # test Laplacian
    du21 = L2d * ue1
    d2err = norm(du21 - ued21, Inf)
    pass = (d2err < 0.75) && (dxerr < 0.1) && (dyerr < 1.e-12)
    pass || println("Discretization test fails")
    return pass
end

"""
lap_test()

Does the FFT invert the discrete Laplacian?
Does the discrete Laplacian pass a simple eigenvalue test.
"""
function lap_test(n, x, z, L2d, fdata)
    randok = rand_test(z, n, L2d, fdata)
    eigok = eig_test(n, x, fdata)
    pass = randok && eigok
    return pass
end


function rand_test(z, n, L2d, fdata)
    n2 = n * n
    z1 = reshape(z, (n2,))
    y1 = L2d * z1
    y = reshape(y1, (n, n))
    mz = fish2d(y, fdata)
    q = reshape(mz, (n2,))
    pass = (norm(q - z1, Inf) < 1.e-12)
    pass || println("rand_test fails, norm =", norm(q - z1))
    return pass
end

function eig_test(n, x, fdata)
    lambda = pi * pi * 5
    efunx = sin.(pi * x)
    efuny = sin.(2 * pi * x)
    efunu = efunx * efuny'
    vfun = fish2d(efunu, fdata)
    pass = (norm(lambda * vfun - efunu, Inf) < 1.e-2)
    pass || println("eig test fails")
    return pass
end
