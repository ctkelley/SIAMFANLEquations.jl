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
nsoloutfd=nsold(heqf!, x0, FS, FPS; pdata = hdata);
nsoloutbos=nsold(heqbos!, x0, FS, FPS; pdata = c);
dbos=norm(nsoloutbos.solution-nsoloutfd.solution)
bosok=dbos<1.e-7
if ~bosok
    println("Bosma and DeRooij test fails in H-equation")
end
nsoloutsp=nsold(heqf!, x0, FS, FPSS, heqJ!; pdata = hdata)
dsp=norm(nsoloutsp.history-nsoloutfd.history)
dsolsp=norm(nsoloutsp.solution-nsoloutfd.solution)
spok = (dsp<1.e-8) && (dsolsp<1.e-6)
if ~spok
    println("Mixed precision test fails in H-equation")
end
#
# change c and use old solution as initial iterate
#
h5=nsoloutfd.solution
setc!(hdata,.7)
nsoldoutfd7=nsold(heqf!, h5, FS, FPS; pdata = hdata)
contok = (nsoldoutfd7.history[4] < 1.e-12)
if ~contok
    println("Update c test fails in H-equation")
end
heqok=spok && bosok && contok
if ~heqok
   println("H-equation Chapter 2 test fails")
end
return heqok
end
