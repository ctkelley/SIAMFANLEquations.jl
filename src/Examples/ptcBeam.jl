"""
ptcBeam(n, maxit, pdt=.01, lambda=20.0; precision=Float64, keepsolhist=false)

Test PTC for systems on the buckling beam problem.
Compare to Newton, which will converge to the unstable solution.
"""
function ptcBeam(
    n,
    maxit,
    pdt = 0.01,
    lambda = 20.0;
    precision = Float64,
    keepsolhist = false,
    jknowsdt = false,
)
    #
    # This is a steady-state computation, so there is no dt in the problem.
    #
    bdata = beaminit(n, 0.0, lambda)
    x = bdata.x
    u0 = x .* (1.0 .- x) .* (2.0 .- x)
    u0 .*= exp.(-10.0 * u0)
    FS = copy(u0)
    FPS = precision.(copy(bdata.D2))
    if jknowsdt
        Jeval = BeamJdt!
    else
        Jeval = BeamJ!
    end
    bout = ptcsol(
        FBeam!,
        u0,
        FS,
        FPS,
        Jeval;
        #        BeamJ!;
        rtol = 1.e-10,
        pdata = bdata,
        pdt0 = pdt,
        maxit = maxit,
        jknowsdt = jknowsdt,
        keepsolhist = keepsolhist,
    )
    if ~jknowsdt
        qout = nsol(FBeam!, u0, FS, FPS, BeamJ!; pdata = bdata, sham = 1)
        return (bout, qout)
    else
        return bout
    end
end

"""
BeamJdt!(FP, FV, U, dt, bdata)
Evaluates the Jacobian + (1/dt) I for PTC. 
"""
function BeamJdt!(FP, FV, U, dt, bdata)
    FP .= BeamJ!(FP, FV, U, bdata)
    FP .= FP + (1.0 / dt) * I
end
