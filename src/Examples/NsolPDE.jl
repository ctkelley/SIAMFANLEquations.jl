"""
NsolPDE(n)

Solve the Elliptic PDE using nsol.jl on an n x n grid. You give me
n and (optionally) sham and resdec and I return the output of nsol.
"""
function NsolPDE(n; sham=1, resdec=.5)
u0=zeros(n*n,)
FV=copy(u0)
pdata=pdeinit(n)
u2d=pdata.uexact
exsol=reshape(u2d,(n*n,))
J=copy(pdata.D2)
hout=nsol(pdeF!, u0, FV, J, pdeJ!; rtol=1.e-7, atol=1.e-10, 
          jfact=lu, pdata=pdata, sham=sham, resdec=resdec)
return hout
end
