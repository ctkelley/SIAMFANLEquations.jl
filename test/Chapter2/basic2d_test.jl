"""
basic2d_test()

Test nsold with the simple 2D problem.

"""
function basic2d_test()
x0 = ones(2,1);
fv = zeros(2,1);
jv = zeros(2,2);
jsv = zeros(Float32,2,2);
# local convergence testing
#
# single vs double Jacobian
#
nout=nsold(basic2d!,x0,fv,jv; rtol=1.e-10);
sout=nsold(basic2d!,x0,fv,jsv);
dss=norm(nout.solution-sout.solution)
hss=norm(nout.history-sout.history)
singleok=(norm(dss) < 1.e-7) && (norm(hss) < 1.e-7)
if ~singleok
    println("Single/Double test fails.")
end
#
# chord vs Newton
#
cout=nsold(basic2d!,x0,fv,jv; solver="chord");
dsc=norm(nout.solution-cout.solution)
lch=length(cout.history)
lnh=length(nout.history)
jevals=sum(cout.stats.ijac)
chordok=(dsc < 1.e-6) && (lch == 17) && (lnh == 5) && (jevals == 1)
if ~chordok
    println("Chord/Newton test fails")
end
#
# Analytic vs finite-difference Jacobian
#
eout=nsold(basic2d!,x0,fv,jv,jbasic2d!)
fdok = (norm(eout.history-nout.history) < 1.e-6 ) && 
       (norm(eout.solution-nout.solution) < 1.e-10 ) 
if ~fdok
    println("FD/Analytic test fails.")
end
#
# Shamanskii
#
s1out=nsold(basic2d!,x0,fv,jv; sham=2, rtol=1.e-10);
dout1=norm(s1out.solution - nout.solution)
jevals1=sum(s1out.stats.ijac)
s2out=nsold(basic2d!,x0,fv,jv; sham=2, rtol=1.e-10, resdec=.5);
jevals2=sum(s2out.stats.ijac)
dout2=norm(s2out.solution - nout.solution)
shamok=(dout1 < 1.e-10) && (dout2 < 1.e-10) && (jevals1==4) && (jevals2==3)
if ~shamok
    println("Shamanskii test fails.")
end
#
# Global convergence
#
x0a=[2, .5];
FS=zeros(2,);
FPS=zeros(2,2);
FPSS=zeros(Float32,2,2);
nouta=nsold(simple!, x0a, FS, FPS; keepsolhist=true);
noutb=nsold(simple!, x0a, FS, FPSS, jsimple!; keepsolhist=true);
iarm=nouta.stats.iarm
armok = (iarm[2]==2)
preok = (norm(noutb.solhist - nouta.solhist,Inf) < 1.e-6)
solok = (norm(noutb.solution - nouta.solution,Inf) < 1.e-10)
globok = armok && preok && solok
if ~globok
    println("Global test fails.")
end
return chordok && singleok && fdok && shamok && globok
end


function basic2d!(FV,x)
FV[1]=x[1]*x[1]  -2.0;
FV[2]=exp(x[1]-1) + x[2]*x[2] - 2.0;
end

function jbasic2d!(FV,JV,x)
JV[1,1]=2*x[1] 
JV[1,2]=0.0
JV[2,1]=exp(x[1]-1) 
JV[2,2]=2*x[2]
end
