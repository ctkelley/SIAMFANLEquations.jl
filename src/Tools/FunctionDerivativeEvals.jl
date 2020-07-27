"""
PrepareDerivative(ItRules,x,xm,fc,fm)
"""
function PrepareDerivative(ItRules,x,xm,fc,fm)
newjac=0
newfun=0
fp=ItRules.fp
f=ItRules.f
dx=ItRules.dx
solver=ItRules.solver
if solver == "secant"
   df = (fc - fm) / (x - xm)
   newfun=newfun+1
else
   df = fpeval_newton(x, f, fc, fp, dx)
   newjac=newjac+1
end
return df
end
