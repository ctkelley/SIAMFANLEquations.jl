"""
nk_pde(n)

Solve the Elliptic PDE using nsol.jl on an n x n grid. 
"""
function nk_pde()
# Get some room for the residual
rtol=1.e-7; atol=1.e-10; n=15;
u0=zeros(n*n,); FV=copy(u0);
# Get the precomputed data from pdeinit
pdata=pdeinit(n);
# Storage for the Jacobian-vector products
JV=zeros(n*n,100);
# Call the solver
hout=nsoli(pdeF!, u0, FV, JV; rtol=rtol, atol=atol,
          pdata=pdata, fixedeta=false, maxit=20);
hout2=nsoli(pdeF!, u0, FV, JV, Jvec2d; rtol=rtol, atol=atol, Pvec=Pvec2d,
          pdata=pdata, fixedeta=false, maxit=20, pside="right");
hout3=nsoli(pdeF!, u0, FV, JV, Jvec2d; rtol=rtol, atol=atol, Pvec=Pvec2d,
          pdata=pdata, fixedeta=true, maxit=20, pside="right");
soldiff= (norm(hout3.solution-hout.solution,Inf)
             +norm(hout3.solution-hout.solution,Inf))
solpass=(soldiff < 1.e-6)
histdiff= (hout.history-hout2.history)./hout.history[1]
histpass= (norm(histdiff,Inf) < .1)
cost1=sum(hout.stats.ijac)
cost2=sum(hout2.stats.ijac)
cost3=sum(hout3.stats.ijac)
costpass = (cost1 > 90) && (cost2 > 35) && (cost3 > cost2)
return costpass
end
