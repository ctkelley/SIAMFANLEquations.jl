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
h=1/(n+1)
maindiag=4*ones(n^2,)/(h*h)
mdiag=Pair(0,maindiag);
sxdiag=-ones(n^2-1,)/(h*h);
for iz=n:n:n^2-1
   sxdiag[iz]=0.0;
end
upxdiag=Pair(1,sxdiag);
lowxdiag=Pair(-1,sxdiag);
sydiag=-ones(n^2-n,)/(h*h);
upydiag=Pair(n,sydiag);
lowydiag=Pair(-n,sydiag);
L2d=spdiagm(lowydiag, lowxdiag, mdiag, upxdiag, upydiag);
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
T=fdata.T
(nx,ny) = size(f)
if nx != ny
   error("need a square grid in fish2d")
end
u=sintv(f,fdata)';
u[:]=T\u[:];
u=isintv(u',fdata)
return u
end


"""
uhat=sintv(u, fdata)

Fast sine transform using FFTW 

Preallocated data and using plan_fft!

Nothing fancy.
"""
function sintv(u, fdata)
utmp=fdata.utmp;
FFF=fdata.FFF;
uhat=fdata.uhat
(nx,ny)=size(u);
utmp.=0.0;
@views utmp[2:nx+1,:].=u;
FFF*utmp;
@views uhat.=-imag.(utmp[2:nx+1,:])
return uhat
end


"""
uhat=isintv(u, fdata)

Fast inverse sine transform using FFTW 

Preallocated data and using plan_fft!

Nothing fancy.
"""
function isintv(u, fdata)
utmp=fdata.utmp
IFFF=fdata.IFFF
(nx,ny)=size(u)
utmp .= 0.0
@views utmp[2:nx+1,:].=u;
IFFF*utmp
@views iuhat=4*imag.(utmp[2:nx+1,:])
return iuhat
end

"""
fishinit(n)

Run plan_fft! and plan_ifft! to set up the solver. Do not mess
with this function.
"""
function fishinit(n)
bsize=(2*n+2,n);
zstore = zeros(bsize) * (1.0 + im)
FFF=plan_fft!(zstore,[1]);
IFFF=plan_ifft!(zstore,[1]);
utmp = zeros(bsize) * (1.0 + im)
uhat=zeros(n,n);
T=newT(n)
fdata=(FFF=FFF, IFFF=IFFF, utmp=utmp, uhat=uhat, T=T)
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
