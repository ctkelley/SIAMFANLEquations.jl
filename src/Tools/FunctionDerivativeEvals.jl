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

function UpdateIteration(xt::T, xm, lambda, d, ItRules) where T <: Real
f=ItRules.f
xt = xm + lambda * d
fc = f(xt)
residc=norm(fc)
return (xt, residc, fc)
end
