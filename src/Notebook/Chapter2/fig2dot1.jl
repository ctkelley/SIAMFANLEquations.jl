"""
fig2dot1(inbook=false)
Draw Fig 2.1 in the print book.
"""
function fig2dot1(inbook=false)
x0a=[2, .5];
x0b=[3, 5];
FS=zeros(2,);
FPS=zeros(2,2);
nouta=nsold(simple!, x0a, FS, FPS; keepsolhist=true);
xa=nouta.solhist[1,:];
ya=nouta.solhist[2,:];
noutb=nsold(simple!, x0b, FS, FPS; keepsolhist=true);
xb=noutb.solhist[1,:];
yb=noutb.solhist[2,:];
plot(xa,ya,"k-*")
plot(xb,yb,"k-",xb,yb,"ko")
contourch2()
if ~inbook
title("Fig 2.1 from print book")
end
end


"""
contourch2()

Makes contour plot in Chapter 2
"""
function contourch2()
yr=-5:.01:5
xr=0:.01:5
xmesh=meshgrid(xr,yr)
x=xmesh.x
y=xmesh.y
la=0:.25:.5;
lc=.6:.4:1;
lb=2:4:50;
levels=[la;lc;lb]
zx=x.*x + y.*y .- 2;
zy=exp.(x.-1) + y.*y .-2
z=sqrt.(zx.*zx + zy.*zy)
contour(x,y,z,levels;colors=["black"])
xlabel(L"x_1")
ylabel(L"x_2")
end


"""
duplicate of matlab meshgrid function
"""
function meshgrid(xin,yin)
nx=length(xin)
ny=length(yin)
xout=zeros(ny,nx)
yout=zeros(ny,nx)
for ix=1:ny
    for jx=1:nx
        xout[ix,jx]=xin[jx]
        yout[ix,jx]=yin[ix]
end
end
return (x=xout, y=yout)
end

