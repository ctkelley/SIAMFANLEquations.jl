function ptcKrylovTestC(n = 63)

maxit=100; delta0 = 0.01; lambda = 20.0;
pout1 = ptciBeam()
bdata = beaminit(n, 0.0, lambda);
x = bdata.x; u0 = x .* (1.0 .- x) .* (2.0 .- x); 
u0 .*= exp.(-10.0 * u0); FS = copy(u0); FPJV=zeros(n,20);
pout = ptcsoli( FBeam!, u0, FS, FPJV; 
         delta0 = delta0, pdata = bdata, eta = 1.e-2, 
         rtol = 1.e-10, maxit = maxit, Pvec = PreCondBeam);
delsol=norm(pout.solution-pout1.solution,Inf)
hpass=(length(pout.history) == 25)
solpass=(delsol < 1.e-9)
ptciok = hpass && solpass
end

function PreCondBeam(v, x, bdata)
    J = bdata.D2
    ptv = J\v
end

