"""
ptciBeam(n=63, delta0=1.e-2, PvecKnowsdelta=true, pside = "right")

Solves the buckling beam problem with ptcsoli. You can play
- left/right preconditioning
- pseudo time step dependent preconditioning
- relationship of delta0 to n (hint, it's not mesh-independent)
"""
function ptciBeam(n = 63, delta0 = 1.e-2, PvecKnowsdelta = true, pside = "right")
    lambda = 20.0
    maxit = 1000
    delta0 = 0.01
    PvecKnowsdelta ? Pvec = ptvbeamdelta : Pvec = ptvbeam
    bdata = beaminit(n, 0.0, lambda)
    x = bdata.x
    u0 = x .* (1.0 .- x) .* (2.0 .- x)
    u0 .*= exp.(-10.0 * u0)
    FS = copy(u0)
    FPJV = zeros(n, 20)
    pout = ptcsoli(
        FBeam!,
        u0,
        FS,
        FPJV;
        delta0 = delta0,
        pdata = bdata,
        eta = 1.e-2,
        rtol = 1.e-10,
        maxit = maxit,
        Pvec = Pvec,
        PvecKnowsdelta = PvecKnowsdelta,
        pside = pside,
    )
    return pout
end



"""
ptvbeamdelta(v, x, bdata)

Precondition buckling beam problem with delta-aware preconditioner.
"""
function ptvbeamdelta(v, x, bdata)
    delta = bdata.deltaval[1]
    J = bdata.D2 + (1.0 / delta) * I
    ptv = J \ v
end

"""
ptvbeamp(v, x, bdata)

Precondition buckling beam problem with inverse of high-order term.
"""
function ptvbeam(v, x, bdata)
    J = bdata.D2
    ptv = J \ v
end
