"""
heq_lu_test()

Does the H-equation module do what it's supposed to?
"""
function heq_lu_test()
n=32;
c=.5;
FS=ones(n);
x0=ones(n);
FPS=ones(n,n);
FPSS=ones(Float32,n,n);
hdata = heqinit(x0, c)
nsoloutfd=nsolheq(x0, FS, FPS, hdata)
nsoloutbos=nsol(heqbos!, x0, FS, FPS; pdata = c);
dbos=norm(nsoloutbos.solution-nsoloutfd.solution)
bosok=dbos<1.e-7
if ~bosok
    println("Bosma and DeRooij test fails in H-equation")
end
nsoloutsp=nsolheq(x0, FS, FPSS, hdata; diff=:exact)
nsoloutdp=nsolheq(x0, FS, FPS, hdata; diff=:exact)
dsp=norm(nsoloutsp.history-nsoloutfd.history)
ddp=norm(nsoloutdp.history-nsoloutfd.history)
dsolsp=norm(nsoloutsp.solution-nsoloutfd.solution)
dsoldp=norm(nsoloutsp.solution-nsoloutdp.solution)
spok = (dsp<1.e-7) && (dsolsp<1.e-9) && (ddp < 1.e-7) && (dsoldp < 1.e-9)
if ~spok
    println("Mixed precision test fails in H-equation")
end
#
# change c and use old solution as initial iterate
#
h5=nsoloutfd.solution
setc!(hdata,.7)
nsoloutfd7=nsolheq(h5, FS, FPS, hdata)
contok = (nsoloutfd7.history[4] < 1.e-12)
if ~contok
    println("Update c test fails in H-equation")
end
heqok=spok && bosok && contok
if ~heqok
   println("H-equation Chapter 2 test fails")
end
return heqok
end
