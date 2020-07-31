"""
PrepareJac!(FS, FPS::Array{T,2}, x, F!, J!, dx, pdata; 
                    fact! = lu!) where T<:Real

Compute the Jacobian and perform the factorization. If know something
about the Jacobian, you can tell me what factorization to use. 

For example, if your Jacobian is spd, fact!=cholesky! would work well.

"""
#function PrepareJac!(FS, FPS::Array{T,2}, x, F!, J!, dx, pdata; 
#                     fact! = lu!) where T<:Real
function PrepareJac!(FS, FPS::Array{T,2}, x, ItRules) where T<:Real
F! =ItRules.F!
J! =ItRules.J!
dx =ItRules.dx
fact! = ItRules.fact!
pdata=ItRules.pdata
EvalJ!(FS, FPS, x, F!, J!, dx, pdata)
fact!(FPS)
end

"""
EvalF!(FS, x, F!, pdata)

This is a wrapper for the function evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function EvalF!(FS, x, F!, q::Nothing)
        F!(FS, x)
end

function EvalF!(FS, x, F!, pdata)
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
        EvalF!(FY, y, F!, pdata)
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
F! = ItRules.F!
pdata = ItRules.pdata
BLAS.axpy!(lambda,step,x)
EvalF!(FS, x, F!, pdata)
resnorm=norm(FS) 
return(x, FS, resnorm)
end
