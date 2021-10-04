"""
ci_pde_aa()
Duplicate part of the data for Figure 4.2 in the book.
"""
function ci_pde_aa()
n=31;
m=10;
pdata=pdeinit(n);
Vstore=zeros(n*n,3*m+3);
VstoreS=zeros(n*n,2*m+4);
aout=PDE_aa(n, m; Vstore=Vstore, pdata=pdata);
aoutS=PDE_aa(n, m; Vstore=VstoreS, pdata=pdata);
# Same results with low storage mode?
alphaS= reldiff(aout.stats.alphanorm,aoutS.stats.alphanorm)
condS= reldiff(aout.stats.condhist,aoutS.stats.condhist)
histS = norm(aoutS.history-aout.history,Inf)
pdeerrS=  condS + alphaS + histS
aout.idid || println("pde solver failed")
(aout.errcode == 0) || println("wrong error code in pde")
(pdeerrS < 1.e-8) || println("different stats   ",condS,"  ",alphaS,"  ",histS)
(length(aout.history)==21) || println("history length wrong")
aa_ok = aout.idid && (aout.errcode==0) && (length(aout.history)==21) && (pdeerrS < 1.e-8)
aa_ok && println("pde succeeds")
return aa_ok
end
