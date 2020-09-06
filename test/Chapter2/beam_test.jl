"""
function beam_test()
Test the time-dependent and steady state beam problem.

"""
function beam_test()
#
dt=.02; n=20; stepnum=5
(t, se, xe, fhist, fhistt) =  beamivp(n, dt, stepnum)
beamtdout= (length(fhist)==6) && (norm(fhistt,Inf) < 5.e-5)
return beamtdout
end

"""
beamivp(n, dt, stepnum)
Solve the time-dependent beam problem. Return the iteration history
for the figure and the tables.

"""
function beamivp(n, dt, stepnum)
#
# Set up the initial data for the temporal integration and
# the figure.
#
    bdata = beaminit(n, dt)
    FB = zeros(n)
    FR = zeros(n)
    zd = zeros(n)
    zr = zeros(n - 1)
    JB = Tridiagonal(zr, zd, zr)
    x = bdata.x
    un = (1.e-4) * x .* (1.0 .- x)
    bdata.UN .= un
    nout = []
    solhist = zeros(n, stepnum + 1)
    solhist[:, 1] .= un
    fhist = []
    fhistt = []
    idid = true
    idt = 1
    fx = FBeam!(FR, un, bdata)
    fxn = norm(fx, Inf)
    fxt = FBeamtd!(FR, un, bdata)
    fxtn = norm(fxt)
    push!(fhist, fxn)
    push!(fhistt, fxtn)
    #
    # Take a stepnum time steps and accumulate the data for the book.
    # The integration will terminate prematurely if the nonlinear solve fails.
    # This can happen if your time step is too large and/or your 
    # predictor is poor.
    #
    # I have tuned the time step to make the solver happy and 
    # we are getting close to steady state.
    #
    while idt <= stepnum && idid && fxn > 1.e-12
        nout = nsold(
            FBeamtd!,
            un,
            FB,
            JB,
            BeamtdJ!;
            pdata = bdata,
            atol = 1.e-6,
            rtol = 1.e-6,
            maxit = 3,
            solver = "chord",
        )
        idid = nout.idid
        un = nout.solution
        solhist[:, idt+1] .= un
        bdata.UN .= un
        idt += 1
        fx = FBeam!(FR, un, bdata)
        fxn = norm(fx, Inf)
        push!(fhist, fxn)
        push!(fhistt, nout.history[end])
    end
    t = dt * collect(1:1:idt)
    zp=zeros(idt,);
    se=[zp solhist[:,1:idt]' zp]'
    xe=[0.0 x' 1.0]'
    return (t=t, se=se, xe=xe, fhist=fhist, fhistt=fhistt)
end


