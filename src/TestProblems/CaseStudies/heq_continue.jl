function heq_continue(n=100;version=3)
    #
    # Original form: heqfv1!. c is the parameter
    # PitCon form: heqfv2!. x[n] is the parameter
    #
    FPS = zeros(n, 20);
    FS = zeros(n);
    initname="solutionv"*string(version)*"_init"
    initarray=[solutionv1_init, solutionv2_init, 
               solutionv3_init]
    (FFUN, fdata, pval, nval, xin, x, x0, xold, xdot, 
     bif_update, setlam,lambda, dlam, lambdamax) = initarray[version](n,FPS,FS);
    qdata = (fdata = fdata, FS = FS, FPS = FPS, dlam = dlam, xix=xin, 
      xold=xold, xdot=xdot, lambdamax = lambdamax, bif_update=bif_update, 
      setlam=setlam);
    (pval, nval, x, lambdaz) = knl_continue(FFUN, qdata, pval, 
                nval, x, x0, lambda)
    return (pval=pval, nval=nval, x=x, lambdaz=lambdaz)
end

function solutionv2_init(n,FPS=[],FS=[])
return nothing
end

function solutionv1_init(n,FPS=[],FS=[])
    x0 = ones(n)
    x = copy(x0)
    xold = copy(x0)
    xdot = copy(x0)
    xin = copy(x0)
    lambda=0.0
    dlam = 0.01
    lambdamax = 1.0
    pval = [lambda]
    nval = [norm(x, 1) / n]
    hdata = heqinit(x0, 0.5)
    FFUN = heqf!
#
    function bif_update_1!(pval, nval, x, lambda)
        c = lambda
        n = length(x)
        push!(nval, norm(x, 1) / n)
        push!(pval, c)
    end
#
    function setlam_v1!(qdata, lambda, xdot=[], xold=[])
    hdata=qdata.fdata
    c=lambda
    setc!(hdata,c)
    end
#
    return (
        FFUN = FFUN,
        fdata = hdata,
        pval = pval,
        nval = nval,
        xin = xin,
        x = x,
        x0 = x0,
        xold= xold,
        xdot = xdot,
        bif_update = bif_update_1!,
        setlam = setlam_v1!,
        lambda = lambda,
        dlam = dlam,
        lambdamax = lambdamax
    )
end

function solutionv3_init(n,FPS=[], FS=[])
# Solution at s=0, which we will not compute
# This is the pseudo-arclength version, so x = (H, c)
pval=[0.0]
nval=[1.0]
FFUN=heqfv3!
#
lambda=0.0
dlam = 0.1
#
#
lambdamax = 150.0
#
# Set up the (bogus) precomputed data
#
x0=ones(n)
znew = ones(n)
xin = ones(n-1)
@views xin .= x0[1:n-1]
hdata=heqinit(xin,dlam)
#
# Now compute an honest solution to start the continuation
# we need at least two points on the path before we can come up 
# with xdot. The plan is to solve the equation and then approximate
# xdot vi xdot = (xc - xold)/ds
#
FST=zeros(n-1); FSTP=zeros(n-1,20)
nout=nsoli(heqf!, xin, FST, FSTP; pdata=hdata)
#
@views znew[1:n-1] .= nout.solution
znew[n]=dlam; push!(pval,dlam);
@views nrm=norm(znew[1:n-1],1)/(n-1)
push!(nval,nrm)
zold=ones(n); zold[n]=0.0
xdot = (znew - zold)/dlam
xold = znew
x0 = 2.0*znew - zold
fdata=(hdata=hdata, xold=xold, xdot = xdot, dlam=dlam, xin=xin)
    return (
        FFUN = FFUN,
        fdata = fdata,
        pval = pval,
        nval = nval,
        xin = xin,
        x = znew,
        x0 = x0,
        xold=xold,
        xdot = xdot,
        bif_update = bif_update_3!,
        setlam = setlam_v3!,
        lambda = lambda,
        dlam = dlam,
        lambdamax = lambdamax
    )


end

    function setlam_v3!(qdata, lambda, xdot, xold)
    qdata.fdata.xdot .= xdot
    qdata.fdata.xold .= xold
    end

    function bif_update_3!(pval, nval, x, lambda)
    n=length(x);
    c=x[n];
    push!(pval,c)
    @views nrm=norm(x[1:n-1],1)/(n-1);
    push!(nval,nrm);
    end


#
# Make the input z = ((x, c) , s)
#
function heqfv3!(F, z, fdata)
np=length(z)
n=np-1
hdata=fdata.hdata
dlam = fdata.dlam
zdot = fdata.xdot
zold = fdata.xold
xin = fdata.xin
cdot = zdot[np]
cold = zold[np]
@views xold = zold[1:n]
@views xdot = zdot[1:n]
@views xin .= z[1:n]
c=z[np];
setc!(hdata,c)
@views FST = F[1:n]
FST = heqf!(FST, xin, hdata)
dx = xin
dx .-= xold
Nval = 100.0*(dot(xdot,dx)/n) + cdot*(c - cold) - dlam
F[np] = Nval
return F
end


