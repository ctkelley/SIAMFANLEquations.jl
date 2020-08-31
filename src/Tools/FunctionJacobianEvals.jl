"""
PrepareJac!(FPS::Array{T,2}, FS, x, ItRules) where T<:Real

Compute the Jacobian and perform the factorization. If know something
about the Jacobian, you can tell me what factorization to use. 

For example, if your Jacobian is spd, fact!=cholesky! would work well.

"""
function PrepareJac!(FPS, FS, x, ItRules) 
F! =ItRules.f
J! =ItRules.fp
dx =ItRules.dx
fact = ItRules.fact
pdata=ItRules.pdata
EvalJ!(FS, FPS, x, F!, J!, dx, pdata)
TF=fact(FPS)
return TF
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

function klfact(A::BandedMatrix)
TF=qr(A)
end

function klfact(A::Tridiagonal)
TF=A
end

function klfact(A)
TF = lu(A)
end


"""
EvalF!(F!, FS, x, pdata)

This is a wrapper for the function evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function EvalF!(F!, FS, x, q::Nothing)
        F!(FS, x)
end

function EvalF!(F!, FS, x, pdata)
        F!(FS, x, pdata)
end


"""
JV!(FS, FPS, x, J!, pdata)

This is a wrapper for the Jacobian evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function JV!(FS, FPS, x, J!, pdata)
        J!(FS, FPS, x, pdata) 
end

function JV!(FS, FPS, x, J!, q::Nothing)
        J!(FS, FPS, x) 
end

"""
EvalJ!(FS, FPS, x, F!, J!, dx, pdata)

evaluates the Jacobian before the factorization in PrepareJac!

"""

function EvalJ!(FS, FPS, x, F!, J!, dx, pdata)
    if J! != diffjac!
        JV!(FS, FPS, x, J!, pdata)
    else
        diffjac!(FS, FPS, F!, x, dx, pdata)
    end
end

"""
   diffjac!(FS, FPS::Array{T,2}, F!, x, dx, pdata) where T <: Real

Computes a finite-difference dense and unstructured Jacobian.
This is not something an user wants to mess with. Look at the 
docstrings to nsold to see more details.


Nothing much to see here. Move along.
"""
function diffjac!(FS, FPS::Array{T,2}, F!, x, dx, pdata) where T <: Real
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

