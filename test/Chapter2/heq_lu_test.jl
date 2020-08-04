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
hdata = heqinit(x0, n, c, Float64)
nsoloutfd=nsold(heqf!, x0, FS, FPS; pdata = hdata);
nsoloutbos=nsold(heqbos!, x0, FS, FPS; pdata = c);
dbos=norm(nsoloutbos.solution-nsoloutfd.solution)
bosok=dbos<1.e-7
nsoloutsp=nsold(heqf!, x0, FS, FPSS, heqJ!; pdata = hdata)
dsp=norm(nsoloutsp.history-nsoloutfd.history)
dsolsp=norm(nsoloutsp.solution-nsoloutfd.solution)
spok = (dsp<1.e-8) && (dsolsp<1.e-6)
heqok=spok && bosok
if ~heqok
   println("H-equation Chapter 2 test fails")
end
return heqok
end
