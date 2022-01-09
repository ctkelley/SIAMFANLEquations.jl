"""
heat_fixed!(theta, thetain, hn_data)

Fixed point map for the conductive-radiative heat transfer problem.
"""
function heat_fixed!(theta, thetain, hn_data)
epsl=1.0; epsr=1.0; thetal=1.0; thetar=0.0;
sn_data=hn_data.sn_data
nx=length(thetain)
theta .= thetain
source = sn_data.tmphf
source .*= 0.0
rhsd2 = hn_data.rhsd2
bcfix=hn_data.bcfix
D2=hn_data.D2
Nc=hn_data.Nc
omega=hn_data.omega
source .= theta
source .^=4
source .*= (1.0-omega)
ltol=1.e-12
flux=flux_solve(source, hn_data, ltol)
@views copy!(rhsd2,flux[2:nx-1])
rhsd2 .*= (1.0 - omega)
@views axpy!(-2.0, source[2:nx-1], rhsd2)
pn=1.0/(2.0*Nc)
rhsd2 .*= pn
ldiv!(D2,rhsd2)
theta[1]=0.0; theta[nx]=0.0;
@views theta[2:nx-1] .= rhsd2
axpy!(1.0, bcfix, theta)
return theta
end

"""
heat_init(nx, na, thetal, thetar, omega, tau, Nc)

Set up the conductive-radiative heat transfer problem
"""
function heat_init(nx, na, thetal, thetar, omega, tau, Nc)
np=nx-2
D2M=Lap1d(np)
D2 = ldlt(D2M)
rhsd2=zeros(np)
h=tau/(nx-1.0)
kl_store=kstore(nx,"gmres")
xv=collect(0:h:tau)
bcfix=thetal .+ (thetar - thetal)*xv
sn_data=sn_init(nx, na, x -> omega, tau, thetal, thetar)
hn_data=(sn_data = sn_data, bcfix = bcfix, D2=D2,
         rhsd2 = rhsd2, omega=omega, Nc=Nc,
         kl_store=kl_store)
return hn_data
end


"""
sn_init(nx, na2, fs, tau, vleft, vright; siewert=false)

I pass a named tuple of precomputed and preallocated data to
all the functions and solvers. 

The input to this is obvious stuff.

nx = number of spatial grid points

na2 = number of angles. The angular mesh is (na2/2) Gaussian quadaratures
      on [-1,0) and (0,1]

fs:function ; scattering coefficient is fs(x)


Boundary conditions for the transport problem are constant vectors
filled with vleft/vright.

phi_left, phi_right = ones(na2/2) * vleft/vright
"""
function sn_init(nx, na2, fs, tau, vleft, vright; siewert=false)
#
# Set up the quadrature rule in angle
#
    if siewert
    #
    # I don't need the weights to make tables, but I need
    # to return something.
    #
    angles = [-.05; collect(-.1:-.1:-1.0); 0.05; collect(0.1:0.1:1.0)]
    weights = angles
    else
    (angles, weights) = hard_gauss()
#    (angles, weights) = sn_angles(na2)
    end
    na = floor(Int, na2 / 2)
#
# scattering coefficient
#
    dx = tau / (nx - 1)
    x = collect(0:dx:tau)
    c = fs.(x)
#
# Preallocated storage for intermediate results
#
    phi0 = zeros(nx)
    tmpf = zeros(nx)
    tmp1 = zeros(nx)
    tmphf = zeros(nx)
    rhsg = zeros(nx)
    ptmp = zeros(na)
#
# Preallocates storage for source iteration
#
    psi_left=vleft*ones(na)
    psi_right=vright*ones(na)
    psi = zeros(na2, nx)
    source_average=zeros(nx-1)
    source_total=zeros(nx)
#
# Preallocated storage for the GMRES solve
#
    V=zeros(nx,13)
#
    return sn_data = (
        c = c,
        dx = dx,
        psi = psi,
        angles = angles,
        weights = weights,
        phi0 = phi0,
        tmp1 = tmp1,
        tmpf = tmpf,
        tmphf = tmphf,
        rhsg=rhsg,
        source_average=source_average,
        source_total=source_total,
        nx = nx,
        ptmp = ptmp,
        psi_left = psi_left,
        psi_right = psi_right,
        V = V
    )
end


function hard_gauss()
#
# Return the weights/nodes for double 20 pt gauss
#
m=40
ri=zeros(40)
wi=zeros(40)
r=zeros(40)
w=zeros(40)
       ri[20]=.993128599185095;
       ri[19]=.963971927277914;
       ri[18]=.912234428251326;
       ri[17]=.839116971822218;
       ri[16]=.746331906460151;
       ri[15]=.636053680726515;
       ri[14]=.510867001950827;
       ri[13]=.373706088715420;
       ri[12]=.227785851141645;
       ri[11]=.076526521133497;
       wi[20]=.017614007139152;
       wi[19]=.040601429800387;
       wi[18]=.062672048334109;
       wi[17]=.083276741576705;
       wi[16]=.101930119817240;
       wi[15]=.118194531961518;
       wi[14]=.131688638449177;
       wi[13]=.142096109318382;
       wi[12]=.149172986472604;
       wi[11]=.152753387130726;
       for i=1:10,
          ri[i]=-ri[21-i];
          wi[i]=wi[21-i];
      end
#      mm=m/2;
      mm=floor(Int, m / 2)
      for i=1:mm
       r[i+mm]=(1.0 + ri[i])*.5;
       w[i+mm]=wi[i]*.5;
#       r[i]=(ri[i] - 1.0)*.5;
       r[i]=-r[i+mm]
       w[i]=wi[i]*.5;
       end
return (r,w)
end


"""
sn_angles(na2=20)

Get double Gauss nodes and weights for SN
This function uses FastGaussQuadrature
"""
function sn_angles(na2 = 20)
    na = floor(Int, na2 / 2)
    2 * na == na2 || error("odd number of angles")
    baseangles, baseweights = gausslegendre(na)
    posweights = baseweights * 0.5
    negweights = copy(posweights)
    posangles = (baseangles .+ 1.0) * 0.5
    negangles = -copy(posangles)
    weights = [negweights; posweights]
    angles = [negangles; posangles]
    angles, weights
end


"""
flux_solve(source, hn_data, tol)

Solve the transport equation with the source from the heat
conduction problem. The output is what kl_gmres returns, so
the solution is kout.sol
"""
function flux_solve(source, hn_data, tol)
sn_data=hn_data.sn_data
b=getrhs(source, sn_data)
kl_store=hn_data.kl_store
kout=kl_gmres(sn_data.phi0,b,AxB,sn_data.V,tol; pdata=sn_data,
              kl_store=kl_store)
return kout.sol
end

function AxB(flux, sn_data)
nx=length(flux)
angles=sn_data.angles
na2=length(angles)
na = floor(Int, na2 / 2)
#tmp1=zeros(nx)
#tmpf=zeros(nx)
tmpf = sn_data.tmpf
tmp1 = sn_data.tmp1
tmp1 .*= 0.0
tmpf.=flux
tmp2=zeros(na)
tmpf=source_iteration!(tmpf,tmp2,tmp2,tmp1,sn_data)
axpy!(-1.0, flux, tmpf)
tmpf .*= -1.0
return tmpf
end

function getrhs(source,sn_data)
nx=sn_data.nx
#rhs=zeros(nx)
rhs=sn_data.rhsg
rhs .*= 0.0
angles=sn_data.angles
na2=length(angles)
na = floor(Int, na2 / 2)
rhs=source_iteration!(rhs, sn_data.psi_left, sn_data.psi_right, source, sn_data)
return rhs
end

function source_iteration!(flux,psi_left,psi_right,source,sn_data)
psi=sn_data.psi
psi=transport_sweep!(psi, flux, psi_left, psi_right, source, sn_data)
weights=sn_data.weights
nx=sn_data.nx
na2=length(weights)
    #
    # Take the 0th moment to get the flux.
    #
    g=reshape(flux,1,nx)
    wt=reshape(weights,1,na2)
    mul!(g, wt, psi)
    return flux
end

"""
transport_sweep!(psi, phi, psi_left, psi_right, source, sn_data)

Take a single transport sweep.
"""
function transport_sweep!(psi, phi, psi_left, psi_right, source, sn_data)
    angles = sn_data.angles
    #
    c = sn_data.c
    dx = sn_data.dx
    #
    na2 = length(angles)
    na = floor(Int, na2 / 2)
    nx = length(phi)
    source_average = sn_data.source_average
    source_total= sn_data.source_total
    copy!(source_total,phi)
    source_total .*= .5
    source_total .*= c
    axpy!(1.0, source, source_total)
    @views copy!(source_average, source_total[2:nx])
    @views source_average .+= source_total[1:nx-1]
    source_average .*= .5
    @views forward_angles = angles[na+1:na2]
    @views backward_angles = angles[1:na]
    vfl = (forward_angles / dx) .+ 0.5
    vfl = 1.0 ./ vfl
    vfr = (forward_angles / dx) .- 0.5
    psi .*= 0.0
    @views psi[1:na, nx] .= psi_right
    @views psi[na+1:na2, 1] .= psi_left
    #
    # Forward sweep
    #
    @views for ix = 2:nx
        copy!(psi[na+1:na2, ix], psi[na+1:na2, ix-1])
        psi[na+1:na2, ix] .*= vfr
        psi[na+1:na2, ix] .+= source_average[ix-1]
        psi[na+1:na2, ix] .*= vfl
    end
    #
    # Backward sweep
    #
    @views for ix = nx-1:-1:1
        copy!(psi[1:na, ix], psi[1:na, ix+1])
        psi[1:na, ix] .*= vfr
        psi[1:na, ix] .+= source_average[ix]
        psi[1:na, ix] .*= vfl
    end
    return psi
end
