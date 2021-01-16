"""
PDE_Tools

This file has the operators I need for the PDE example. They
live in a separate file to make the CI easier for me to organize.
"""
# Famous sparse matrices
"""
Dx2d(n)

returns x partial on n x n grid.
Unit square, homogeneous Dirichlet BC
"""
function Dx2d(n)
h=1/(n+1);
ssdiag=ones(n^2-1,)/(2*h);
for iz=n:n:n^2-1
ssdiag[iz]=0.0;
end
updiag=Pair(1,ssdiag);
lowdiag=Pair(-1,-ssdiag);
Dx=spdiagm(lowdiag,updiag);
return Dx
end

"""
Dy2d(n)

returns y partial on n x n grid.
Unit square, homogeneous Dirichlet BC
"""
function Dy2d(n)
h=1/(n+1);
ssdiag=ones(n^2-n,)/(2*h);
updiag=Pair(n,ssdiag);
lowdiag=Pair(-n,-ssdiag);
Dy=spdiagm(lowdiag,updiag);
return Dy
end

"""
Lap2d(n)

returns the negative Laplacian in two space dimensions
on n x n grid.

Unit square, homogeneous Dirichlet BC
"""
function Lap2d(n)
# hm2=1/h^2
hm2=(n+1.0)^2;
maindiag=fill(4*hm2,(n^2,));
sxdiag=fill(-hm2,(n^2-1,));
sydiag=fill(-hm2,(n^2-n,));
for iz=n:n:n^2-1
   sxdiag[iz]=0.0;
end
L2d=spdiagm(-n => sydiag, -1 => sxdiag, 0=> maindiag, 
             1 => sxdiag, n => sydiag);
return L2d
end


"""
u=fish2d(f, fdata)

Fast Poisson solver in two space dimensions.
Same as the Matlab code.
Unit squre + homogeneous Dirichlet BCs.

Grid is nx by nx

You give me f as a two-dimensional vector f(x,y).
I return the solution u.
"""
function fish2d(f, fdata)
u=fdata.utmp
v=fdata.uhat
T=fdata.T
ST=fdata.ST
(nx,ny) = size(f)
nx == ny || error("need a square grid in fish2d")
u.=f; u=ST*u; u=u'
u1=reshape(u,(nx*nx,))
u1 .= T\u1;
u = reshape(u1,(nx,nx))
#u .= u'; u = ST*u; u ./=(2*nx+2);
v .= u'; v = ST*v; v ./=(2*nx+2);
#return u
return v
end


"""
fishinit(n)

Run FFTW.plan_r2r to set up the solver. Do not mess
with this function.
"""
function fishinit(n)
#
# Get the sine transform from FFTW. This is faster/better/cleaner
# than what I did in the Matlab codes.
#
zstore = zeros(n,n)
ST = FFTW.plan_r2r!(zstore, FFTW.RODFT00, 1);
uhat=zeros(n,n);
T=newT(n)
fdata=(ST=ST, uhat=uhat, utmp=zstore, T=T)
return fdata
end

"""
T = newT(n)

Builds the n^2 x n^2 sparse tridiagonal matrix for
the 2D fast Poisson solver.
"""
function newT(n)
N=n*n;
h=1/(n+1);
x=h:h:1-h;
h2=1/(h*h);
LE=2*(2 .- cos.(pi*x))*h2;
fn=ones(N-1,)*h2; gn=ones(N-1,)*h2; dx=zeros(N,);
for k=1:n-1
    fn[k*n]=0.0;
    gn[k*n]=0.0;
    dx[(k-1)*n+1:n*k]=LE[k]*ones(n,);
end
dx[(n-1)*n+1:n*n]=LE[n]*ones(n,);
T=Tridiagonal(-fn, dx, -gn);
return T
end

"""
Use fish2d and reshape for preconditioning.
"""
function Pfish2d(v,fdata)
u2=fdata.uhat
n2=length(v)
n=Int(sqrt(n2))
(n*n == n2) || error("input to Pfish2d not a square array")
v2=reshape(v,(n,n))
u2.=fish2d(v2,fdata)
u=reshape(u2,(n2,))
return u
end

