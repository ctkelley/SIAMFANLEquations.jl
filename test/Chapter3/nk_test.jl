"""
nk_test()

CI for nsoli
"""
function nk_test()
passsimple=nksimple()
passsimple || println("nksimple fails")
nkpass = passsimple
return nkpass
end
"""
nksimple()

Test nsoli with the simple 2D problem and line search failure and success.
"""
function nksimple()
x0=[2.0;.5]
FPS=zeros(2,4)
FPJ=zeros(2,2)
FS=copy(x0)
#
# For the easy problem we will do analytic Jacobians for
# Newton and forward difference directional derivatives for Newton-GMRES
#
dout=nsol(simple!, x0, FS, FPJ, jsimple!; sham=1, keepsolhist=true)
kout=nsoli(simple!, x0, FS, FPS ;eta=1.e-10, lmaxit=2, keepsolhist=true)
dsolhist=norm(kout.solhist-dout.solhist,Inf)
shpass=(dsolhist < 1.e-7)
shpass || println("solhist compare fails in nksimple")
vconverge = krstest(dout,kout)
#
# For the stagnating problem we will do analytic Jacobians for
# Newton and analytic Jacobian-vector products for Newton-GMRES
#
x0=[3.0;5.0]
dout=nsol(simple!, x0, FS, FPJ, jsimple!; sham=1)
kout=nsoli(simple!, x0, FS, FPS, Jvsimple;eta=1.e-10, lmaxit=2)
vdiverge = krstest(dout,kout)
return vconverge && vdiverge && shpass
end

function krstest(dout,kout)
hdiff=norm(kout.history-dout.history,Inf)
hpass = (hdiff < 5.e-7)
hpass || println("history compare fails in nksimple")
#
adiff=kout.stats.iarm-dout.stats.iarm
apass = (sum(adiff)==0)
apass || println("line search compare fails in nksimple")
#
fdiff=kout.stats.ifun-dout.stats.ifun
fpass = (sum(fdiff)==0)
fpass || println("function value compare fails in nksimple")
#
soldiff = kout.solution - dout.solution
solpass = (norm(soldiff,Inf) < 1.e-11)
solpass || println("solution compare fails in nksimple")
krpass = (fpass && apass && hpass && solpass)
end

function Jvsimple(v, FS, x)
#println(typeof(v))
JM=zeros(2,2)
jsimple!(JM,FS,x)
Jvec=JM*v
return Jvec
end
