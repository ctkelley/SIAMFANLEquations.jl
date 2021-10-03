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

"""
PDE_aa(n=31, m=3; Vstore=Float64[], pdata=nothing)
Solve preconditioned convection-diffusion equation with hardwired left
preconditioner using Anderson acceleration.

If you're putting this in a loop, you should allocate Vstore to
zeros(n*n,2*(mmax+1)). Otherwise I will make a decision for you and
allocate for Vstore with each call to this function. The story on
pdata is the same. If you are calling this several times with the
same value of n, build pdata outside the call.
"""
function PDE_aa(n=31, m=3; Vstore=Float64[], pdata=nothing)
#
# Process Vstore and pdata
#
(pdata != nothing) || (pdata=pdeinit(n));
(length(Vstore)>0) || (Vstore=zeros(n*n,3*m+2));
(mv, nv)  = size(Vstore)
dimvtest= ((mv == n*n) && (nv >= 2*(m + 1)))
dimvtest || error("Vstore too small")
#
# Call aasol and return the results.
#
u0=zeros(n*n);
rtol = 1.e-8;
atol = 1.e-8;
aout = aasol(hardleftFix!, u0, m, Vstore; pdata=pdata, maxit=40, rtol=rtol,
             atol=atol);
return aout
end

