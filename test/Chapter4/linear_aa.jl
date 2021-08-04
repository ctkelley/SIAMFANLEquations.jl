"""
linear_aa()

Test aasol.jl for a two-dimensional linear problem
"""
function linear_aa()
maxit = 10
maxm = 2
vdim = 2 * (maxm + 1)
Vstore = zeros(2, vdim)
eigs=[.1, .5]
xstar=ones(2,)
pdata=makeLinpdata(eigs)
m=0
#
# Test for termination on entry
#
x0=[1.0, 1.0]; m=0;
aout = aasol(GLin!, x0, m, Vstore;
     rtol = 1.e-10, pdata=pdata, maxit=maxit, keepsolhist=true);
tflag = (aout.errcode===-1) && aout.idid
tflag || println("Failure in aasol terminate on entry test.")
#
# Test for failure to converge
#
x0=[2.0, 10.0]
m=0
maxit=10
aout = aasol(GLin!, x0, m, Vstore;
     rtol = 1.e-10, pdata=pdata, maxit=maxit, keepsolhist=true);
failflag = ~aout.idid && (aout.errcode == 10)
failflag || println("Linear iteration failure in aasol test fails.")
#
# Test for convergence is two iterations.
#
m=2
aout = aasol(GLin!, x0, m, Vstore;
     rtol = 1.e-10, pdata=pdata, maxit=maxit, keepsolhist=true);
termflag=(length(aout.history) == 4) && 
         (norm(aout.solution - xstar,Inf) < 1.e-14)
termflag || println("Terminate in two aa iterations test fails.")
#
# Now set the eigenvalues to [2.0, 10.0] and beta=-1/9
#
eigs=[2.0, 10.0]
pdata=makeLinpdata(eigs)
beta=-1.0/9.0
maxit=10
m=1
aout1 = aasol(GLin!, x0, m, Vstore;
     rtol = 1.e-10, pdata=pdata, maxit=maxit)
aout2 = aasol(GLin!, x0, m, Vstore;
     rtol = 1.e-10, pdata=pdata, maxit=maxit, beta=beta)
bflag= ~aout1.idid && aout2.idid && (length(aout2.history)==8)
return tflag && failflag && termflag && bflag

end

function GLin!(gout,xin,pdata)
M=pdata.M
b=pdata.b
gout=M*xin+b
return gout
end

function GLinBeta!(gout,xin,pdata)
M=pdata.M
b=pdata.b
beta=pdata.beta
gout = Glin!(gout,xin,pdata)
gout .*=beta
gout .+= (1.0-beta)*xin
return gout
end

function makeLinpdata(eigs, beta=1.0)
U=[1 -1; 1 1]./sqrt(2.0)
V=[3 -4; 4 3]./5.0
S=diagm(eigs)
M=U*S*V'
b=(I-M)*ones(2,);
return (M=M, b=b, beta=beta)
end
