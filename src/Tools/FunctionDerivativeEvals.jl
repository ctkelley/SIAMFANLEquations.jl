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


"""
fpeval_newton

Evaluates f' by differences or the user's code.

"""
function fpeval_newton(x, f, fc, fp, h)
    if fp == difffp
        df = difffp(x, f, fc, h)
    else
        df = fp(x) 
    end
    return df 
end 


"""
difffp

forward differencing for scalar equations
"""
function difffp(x, f, fc, h)
    df = (f(x + h) - fc) / h
    return df
end
