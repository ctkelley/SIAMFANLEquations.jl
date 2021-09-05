"""
pde_aa()
Duplicate part of the data for Figure 4.2 in the book.
"""
function pde_aa()
n=31;
m=10;
pdata=pdeinit(n);
Vstore=zeros(n*n,2*(m+1));
aout=PDE_aa(n, m; Vstore=Vstore, pdata=pdata);
aa_ok = aout.idid && (aout.errcode==0) && (length(aout.history)==21)
return aa_ok
end
