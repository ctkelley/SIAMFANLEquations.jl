"""
ptctest(n, maxit, dt=.01, lambda=20.0; precision=Float64)

Test PTC for systems on the buckling beam problem.
"""
function ptctest(n, maxit, dt=.01, lambda=20.0; precision=Float64)
bdata=beaminit(n,dt,lambda)
x=bdata.x
u0=x.*(1.0 .- x).^2
u0=x.*(1.0 .- x).*(2.0 .- x)
peak=exp.(-10.0*u0)
u0 .*= peak
FS=copy(u0)
FPS=precision.(copy(bdata.D2))
bout=ptcsol(FBeam!, u0, FS, FPS, BeamJ!; 
             rtol=1.e-10, pdata=bdata, dt0=dt, maxit=maxit);
qout=nsold(FBeam!, u0, FS, FPS, BeamJ!; pdata=bdata);
return (bout, qout)
end
