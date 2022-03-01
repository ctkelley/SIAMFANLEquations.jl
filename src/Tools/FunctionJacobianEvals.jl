#
# The functions in this file manage Jacobian evaluations and
# factorizations and function evaluations. The function evaluation bits
# are used in the Newton-Krylov solvers too.
#

"""
For nsoli I use
PrepareJac!(FPS, FS, x, ItRules)
and for ptcsoli
PrepareJac!(FPS, FS, x, ItRules, dt) 

Compute the Jacobian and perform the factorization. If know something
about the Jacobian, you can tell me what factorization to use. 

For example, if your Jacobian is spd, fact=cholesky! would work well.

"""
function PrepareJac!(FPS, FS, x, ItRules)
    F! = ItRules.f
    J! = ItRules.fp
    dx = ItRules.dx
    fact = ItRules.fact
    pdata = ItRules.pdata
    EvalJ!(FPS, FS, x, F!, J!, dx, pdata)
    TF = fact(FPS)
    return TF
end

function PrepareJac!(FPS, FS, x, ItRules, dt)
    dt > 0 || error("dt must be > 0 in PTC")
    F! = ItRules.f
    J! = ItRules.fp
    dx = ItRules.dx
    fact = ItRules.fact
    jknowsdt = ItRules.jknowsdt
    pdata = ItRules.pdata
    EvalJ!(FPS, FS, x, F!, J!, dx, dt, pdata, jknowsdt)
    TF = fact(FPS)
    return TF
end


"""
PrepareJac!(fc, fm::Real, x, xm, ItRules, dt=0)
Scalar equations
"""
function PrepareJac!(fps::Real, fc, x, ItRules, dt = 0)
    newjac = 0
    newfun = 0
    fp = ItRules.fp
    f = ItRules.f
    dx = ItRules.dx
    pdata = ItRules.pdata
    solver = ItRules.solver
    df = fpeval_newton(x, f, fc, fp, dx, pdata)
    dt == 0 || (df += 1.0 / dt)
    newjac = newjac + 1
    return df
end


"""
klfact(A)

Returns the default choice for the factorization unless you tell
me to do something else. QR is the default choice for banded because
that works properly with Float32.

"""
function klfact(A::Array{T,2}) where {T<:Real}
    TF = lu!(A)
end

# The default for sparse is lu. lu! for sparse matrices is 
# too complicated to put in here. You can use lu! if you
# set fact = nofact and manage the factorization in your Jacobian
# evaluation code. You'll also get to manage the storage. There's
# a project in chapter 2 about that.
#
function klfact(A::SparseMatrixCSC{Float64,Int64})
    TF = lu(A)
end


# The default for banded matrices is qr, because I do not trust 
# you to allocate the extra two upper bands so I can use qr!.
# I'm using qr! in the example in Chapter 2. Look at the source
# to see how I did that.
#
function klfact(A::BandedMatrix)
    TF = qr(A)
end

# Default: do nothing.
function klfact(A)
    TF = nofact(A)
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
    FS = F!(FS, x)
    return FS
end

function EvalF!(F!, FS, x, pdata)
    FS = F!(FS, x, pdata)
    return FS
end

function EvalF!(F!, FS::Real, x::Real, q::Nothing)
    FS = F!(x)
    return FS
end

function EvalF!(F!, FS::Real, x::Real, pdata)
    FS = F!(x, pdata)
    return FS
end


"""
If you let me handle dt in PTC 
JV!(FPS, FS, x, J!, pdata)


If you put the (1/dt) * I in the Jacobian yourself
JV!(FPS, FS, x, dt, J!, pdata)

This is a wrapper for the Jacobian evaluation that figures out if
you are using precomputed data or not. No reason to get excited
about this.
"""
function JV!(FPS, FS, x, J!, pdata)
    J!(FPS, FS, x, pdata)
end

function JV!(FPS, FS, x, dt, J!, pdata)
    J!(FPS, FS, x, dt, pdata)
end

function JV!(FPS, FS, x, dt, J!, q::Nothing)
    J!(FPS, FS, x, dt)
end

function JV!(FPS, FS, x, J!, q::Nothing)
    J!(FPS, FS, x)
end


"""
for Newton
EvalJ!(FPS, FS, x, F!, J!, dx, pdata) 

for PTC
EvalJ!(FPS, FS, x, F!, J!, dx, dt, pdata) 

evaluates the Jacobian before the factorization in PrepareJac!

"""

function EvalJ!(FPS, FS, x, F!, J!, dx, dt, pdata, jknowsdt)
    #    if J! != diffjac!
    #        JV!(FPS, FS, x, J!, pdata)
    #    else
    #        diffjac!(FPS, FS, F!, x, dx, pdata)
    #    end
    if jknowsdt
        FPS = JV!(FPS, FS, x, dt, J!, pdata)
    else
        EvalJ!(FPS, FS, x, F!, J!, dx, pdata)
        FPS .= FPS + (1.0 / dt) * I
    end
    return FPS
end

function EvalJ!(FPS, FS, x, F!, J!, dx, pdata)
    if J! != diffjac!
        JV!(FPS, FS, x, J!, pdata)
    else
        diffjac!(FPS, FS, F!, x, dx, pdata)
    end
    return FPS
end


"""
   diffjac!(FPS::Array{T,2}, FS, F!, x, dx, pdata) where T <: Real

Computes a finite-difference dense and unstructured Jacobian.
This is not something an user wants to mess with. Look at the 
docstrings to nsold to see more details.


Nothing much to see here. Move along.
"""
#function diffjac!(FPS::Array{T,2}, FS, F!, x, dx, pdata) where {T<:Real}
function diffjac!(FPS, FS, F!, x, dx, pdata)
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
    return FPS
end

"""
UpdateIteration

Take a trial step. Evaluate the function and the residual norm.
"""
function UpdateIteration(xt::Array{T}, x, FS, lambda, step, ItRules) where {T<:Real}
    F! = ItRules.f
    pdata = ItRules.pdata
    copy!(xt, x)
    BLAS.axpy!(lambda, step, xt)
    EvalF!(F!, FS, xt, pdata)
    resnorm = norm(FS)
    return (xt, FS, resnorm)
end


function UpdateIteration(xt::T, xm, ft, lambda, d, ItRules) where {T<:Real}
    f = ItRules.f
    pdata = ItRules.pdata
    xt = xm + lambda * d
    fc = 0.0
    fc = EvalF!(f, fc, xt, pdata)
    #fc = f(xt)
    residc = norm(fc)
    return (xt, fc, residc)
end


"""
fpeval_newton

Evaluates f' by differences or the user's code.

"""
function fpeval_newton(x, f, fc, fp, h, pdata)
    fps = string(fp)
    df = 0.0
    dps = string(difffp)
    if fps == dps
        df = difffp(x, f, fc, h, pdata)
    else
        df = EvalF!(fp, df, x, pdata)
    end
    return df
end


"""
difffp

forward differencing for scalar equations
"""
function difffp(x, f, fc, h, pdata)
    fph = 0.0
    fph = EvalF!(f, fph, x + h, pdata)
    df = (fph - fc) / h
    #    df = (f(x + h) - fc) / h
    return df
end


"""
test_evaljac(ItRules, itc, newiarm, residratio)

Figures out if it's time to reevaluate and refacto the Jacbian in
Newton's method.
"""
function test_evaljac(ItRules, itc, newiarm, residratio)
    solver = ItRules.solver
    sham = ItRules.sham
    resdec = ItRules.resdec
    evaljacit = (itc % sham == 0 || newiarm > 0 || residratio > resdec)
    chordinit = (solver == "chord") && itc == 0
    evaljac = (evaljacit && solver == "newton") || chordinit || solver == "secant"
end
