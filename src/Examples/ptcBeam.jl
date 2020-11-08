"""
ptcBeam(n, maxit, dt=.01, lambda=20.0; precision=Float64)

Test PTC for systems on the buckling beam problem.
Compare to Newton, which will converge to the unstable solution.
"""
function ptcBeam(n, maxit, dt = 0.01, lambda = 20.0; precision = Float64)
    bdata = beaminit(n, dt, lambda)
    x = bdata.x
    u0 = x .* (1.0 .- x) .* (2.0 .- x)
    u0 .*= exp.(-10.0 * u0)
    FS = copy(u0)
    FPS = precision.(copy(bdata.D2))
    bout = ptcsol(
        FBeam!,
        u0,
        FS,
        FPS,
        BeamJ!;
        rtol = 1.e-10,
        pdata = bdata,
        dt0 = dt,
        maxit = maxit,
    )
    qout = nsol(FBeam!, u0, FS, FPS, BeamJ!; pdata = bdata, sham = 1)
    return (bout, qout)
end
