"""
NsolPDE(n; sham=1, resdec=.5, rtol=1.e-7, atol=1.e-10)

Solve the Elliptic PDE using nsol.jl on an n x n grid. You give me
n and (optionally) sham and resdec and I return the output of nsol.
"""
function NsolPDE(n; sham=1, resdec=.5, rtol=1.e-7, atol=1.e-10)
# Get some room for the residual
u0=zeros(n*n,)
FV=copy(u0)
# Get the precomputed data from pdeinit
pdata=pdeinit(n)
# Storage for the Jacobian
J=copy(pdata.D2)
# Call the solver
hout=nsol(pdeF!, u0, FV, J, pdeJ!; rtol=rtol, atol=atol,
          pdata=pdata, sham=sham, resdec=resdec)
return hout
end
