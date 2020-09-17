"""
PrepareJac!(FPS, FS, x, ItRules, dt=0) 

Compute the Jacobian and perform the factorization. If know something
about the Jacobian, you can tell me what factorization to use. 

For example, if your Jacobian is spd, fact!=cholesky! would work well.

"""
function PrepareJac!(FPS, FS, x, ItRules,dt=0) 
F! =ItRules.f
J! =ItRules.fp
dx =ItRules.dx
fact = ItRules.fact
pdata=ItRules.pdata
EvalJ!(FPS, FS, x, F!, J!, dx, pdata,dt)
TF=fact(FPS)
return TF
end

"""
PrepareJac!(x::Real,xm,fc,fm,ItRules,dt=0)
Scalar equations
"""
function PrepareJac!(x::Real,xm,fc,fm,ItRules,dt=0)
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


"""
klfact(A)

Returns the default choice for the factorization unless you tell
me to do something else. QR is the default choice for banded because
that works properly with Float32.

"""
function klfact(A::Array{T,2}) where T<:Real
TF = lu!(A)
end

# The default is qr, because I do not trust you to allocate
# the extra two upper bands so I can use qr!.
function klfact(A::BandedMatrix)
TF=qr(A)
end

# Default: do nothing.
function klfact(A)
TF=nofact(A)
end

function nofact(A)
TF = A
end


"""
EvalF!(F!, FS, x, pdata)

This is a wrapper for the function evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function EvalF!(F!, FS, x, q::Nothing)
        F!(FS, x)
        return FS
end

function EvalF!(F!, FS, x, pdata)
        F!(FS, x, pdata)
        return FS
end

function EvalF!(F!, FS, x::Real)
         return FS = F!(x)
end


"""
JV!(FPS, FS, x, J!, pdata)

This is a wrapper for the Jacobian evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function JV!(FPS, FS, x, J!, pdata)
        J!(FPS, FS, x, pdata) 
end

function JV!(FPS, FS, x, J!, q::Nothing)
        J!(FPS, FS, x) 
end

"""
EvalJ!(FPS, FS, x, F!, J!, dx, pdata,dt=0)

evaluates the Jacobian before the factorization in PrepareJac!

"""

function EvalJ!(FPS, FS, x, F!, J!, dx, pdata, dt=0)
    if J! != diffjac!
        JV!(FPS, FS, x, J!, pdata)
    else
        diffjac!(FPS, FS, F!, x, dx, pdata)
    end
    if dt>0
       FPS .= FPS + (1.0/dt)*I
    end
end

"""
   diffjac!(FPS::Array{T,2}, FS, F!, x, dx, pdata) where T <: Real

Computes a finite-difference dense and unstructured Jacobian.
This is not something an user wants to mess with. Look at the 
docstrings to nsold to see more details.


Nothing much to see here. Move along.
"""
function diffjac!(FPS::Array{T,2}, FS, F!, x, dx, pdata) where T <: Real
    h = dx * norm(x, Inf) + 1.e-8
    n = length(x)
    y = ones(size(x))
    FY = ones(size(x))
    for ic = 1:n
        y .= x
        y[ic] = y[ic] + h
        EvalF!(F!, FY, y, pdata)
        for ir = 1:n
            FPS[ir, ic] = (FY[ir] - FS[ir]) / h
        end
    end
end

"""
UpdateIteration

Take a trial step. Evaluate the function and the residual norm.
"""
function UpdateIteration(xt::Array{T}, x, FS, lambda, step, ItRules) where T<:Real
F! = ItRules.f
pdata = ItRules.pdata
copy!(xt,x)
BLAS.axpy!(lambda,step,xt)
EvalF!(F!, FS, xt, pdata)
resnorm=norm(FS) 
return(xt, FS, resnorm)
end


"""
PrepareJacobian!(ItRules,x,xm,fc,fm)
Scalar equations
"""
function PrepareJacobian!(ItRules,x,xm,fc,fm)
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

function UpdateIteration(xt::T, xm, ft, lambda, d, ItRules) where T <: Real
f=ItRules.f
xt = xm + lambda * d
fc = f(xt)
residc=norm(fc)
return (xt, fc, residc)
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

