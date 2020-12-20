"""
nsolpde_test(n)

Test elliptic pde with nsol. Newton and Shamanskii. Query
convergence history, accuracy, agreement.
"""
function nsolpde_test(n)
h=1/(n+1); x=collect(h:h:1.0-h); uexact=solexact(x);
ue=reshape(uexact,(n*n,))
houtn=NsolPDE(n)
histn=houtn.history
npass = (length(histn)==7) && (sum(houtn.stats.iarm)==2) &&
          (histn[7]/histn[1] < 1.e-13)
houts=NsolPDE(n; sham=Inf, resdec=.1)
hists=houts.history
spass = (length(hists)==8) && (sum(houts.stats.iarm)==2) &&
          (hists[8]/hists[1] < 1.e-7)
errn=norm(houtn.solution-ue,Inf)
errs=norm(houts.solution-ue,Inf)
delsol= norm(houts.solution-houtn.solution,Inf)
accpass= (errn < 1.e-3) && (errs < 1.e-3) && (delsol < 1.e-8)
return npass && accpass
end

