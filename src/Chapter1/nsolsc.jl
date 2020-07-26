"""
nsolsc(f,x, fp=difffp; rtol=1.e-6, atol=1.e-12, maxit=10,
        solver=:newton, sham=1, armmax=10, resdec=.1,
        armfix=false, printerr=true, keepsolhist=true)

Newton's method for scalar equations. Has most of the features a
code for systems of equations needs.

Input:\n 
f: function\n 
x: initial iterate\n
fp: derivative. If your derivative function is fp, you give me
its name. For example fp=foobar tells me that foobar is your
function for the derivative. The default is a forward difference
Jacobian that I provide.\n


Keyword Arguments (kwargs):\n
rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations\n

solver:\n
Your choices are :newton(default), "secant", or "chord". However, 
you have sham at your disposal only if you chose newton. "chord"
will keep using the initial derivative until the iterate converges,
uses the iteration budget, or the line search fails. It is not the
same as sham=Inf, which is smarter.\n

If you use secant and your initial iterate is poor, you have made
a mistake. I will help you by driving the line search with a finite
difference derivative.\n

sham:\n
This is the Shamanskii method. If sham=1, you have Newton.
The iteration updates the derivative every sham iterations.
The convergence rate has local q-order sham+1 if you only count
iterations where you update the derivative. You need not
provide your own derivative function to use this option. sham=Inf
is chord only if chord is converging well.\n

armmax: upper bound on stepsize reductions in linesearch

resdec: target value for residual reduction. \n
The default value is .1. In the old MATLAB codes it was .5.
I only turn Shamanskii on if the residuals are decreasing
rapidly, at least a factor of resdec, and the line search is quiescent.
If you want to eliminate resdec from the method ( you don't ) then set
resdec = 1.0 and you will never hear from it again.  

armfix:\n
The default is a parabolic line search (ie false). Set to true and
the stepsize will be fixed at .5. Don't do this unless you are doing
experiments for research.

printerr:\n
I print a helpful message when the solver fails. To supress that
message set printerr to false.

keepsolhist:\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

Output:\n
A tuple (solution, functionval, history, stats, idid, solhist) where
history is the vector of residual norms (|f(x)|) for the iteration
and stats is a tuple of the history of (ifun, ijac, iarm), the number
of functions/derivatives/steplength reductions at each iteration.

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian evaluation. I do count them for
the secant method model.

idid=true if the iteration succeeded and false if not.

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true

# Examples
```jldoctest
julia> nsolout=nsolsc(atan,1.0;maxit=5,atol=1.e-12,rtol=1.e-12);

julia> nsolout.history
6-element Array{Float64,1}:
 7.85398e-01
 5.18669e-01
 1.16332e-01
 1.06102e-03
 7.96200e-10
 2.79173e-24
```

```jldoctest
julia> fs(x)=x^2-4.0; fsp(x)=2x;

julia> nsolout=nsolsc(fs,1.0,fsp; maxit=5,atol=1.e-9,rtol=1.e-9);

julia> [nsolout.solhist.-2 nsolout.history]
6×2 Array{Float64,2}:
 -1.00000e+00  3.00000e+00
  5.00000e-01  2.25000e+00
  5.00000e-02  2.02500e-01
  6.09756e-04  2.43940e-03
  9.29223e-08  3.71689e-07
  2.22045e-15  8.88178e-15

```

"""
function nsolsc(
    f,
    x,
    fp=difffp;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    solver = :newton,
    sham = 1,
    armmax = 5,
    resdec = .1,
    armfix = false,
    printerr = true,
    keepsolhist = true
)
    itc = 0
    idid = true
    iline = true
    h = 1.e-7
    #=
     If you like the secant or sham=large methods, I will do a 
     difference Jacobian anyhow if the line search kicks in. 
     You will thank me for this.
     Even if you don't thank me, I will do it anyhow.
    
     If you insist, solver=chord will ignore poor convergence and let
     things go south with no interference. Please don't do that as
     standard procedure and, if you do, don't blame me.
    =#
    if solver == "secant"
        xm = x * 1.0001
        if xm == 0
            xm = 0.0001
        end
        fm = f(xm)
        sham = 1
        fp = difffp
    end
    fc = f(x)
    derivative_is_old = false
    resid = abs(fc)
    newiarm=-1
    iarm=[0]
    ifun=[1]
    ijac=[0]
    ithist = [abs(fc)] 
    if keepsolhist
        solhist = [x]
    end
    tol = rtol * resid + atol
    residratio = 1
    df=0.0
    armstop=true
    while (resid > tol) && (itc < maxit) && armstop
        newjac=0
        newfun=0
        if solver == "secant"
            df = (fc - fm) / (x - xm)
            newfun=newfun+1
        elseif solver == "chord"
            if itc==0
                df = fpeval_newton(x, f, fc, fp, h)
                newjac=newjac+1
            end
            derivative_is_old = false
        else
         if itc % sham == 0 || newiarm > 0 || residratio > resdec 
                df = fpeval_newton(x, f, fc, fp, h)
                newjac=newjac+1
                derivative_is_old = false
            else
                derivative_is_old = true
            end
        end
        xm = x
        fm = fc
        d = -fc / df
        if true == false
        else
        AOUT = armijosc(fc, d, xm, fm, f, h, fp, armmax, 
                        armfix, derivative_is_old)
        armstop=true
#
# If the line search fails and the derivative is current, stop the iteration.
#
        if AOUT.idid == false
            iline = false
            armstop = derivative_is_old
        end
        newjac = newjac + AOUT.newjac
        fc = AOUT.afc
        x = AOUT.ax
        newiarm = AOUT.aiarm
        newfun=newfun+newiarm+1
        derivative_is_old = AOUT.adfo
        d = AOUT.ad
        resid = abs(fc)
        residratio = abs(fc) / abs(fm)
        itc = itc + 1
        newhist = abs(fc)
        end
        if keepsolhist
            newsol = x
            append!(solhist,newsol)
        end
        append!(iarm,newiarm)
        append!(ifun,newfun)
        append!(ijac,newjac)
        append!(ithist,newhist)
    end
    solution = x
    fval = fc
    resnorm = abs(fval)
    resfail = (resnorm > tol) 
    if resfail || iline == false
        idid = false
        resnorm = norm(fval)
    end
    if idid == false && printerr
        NewtonError(resfail, iline, resnorm, maxit, armmax)    
    end
    stats = (ifun=ifun, ijac=ijac, iarm=iarm)
    if keepsolhist
        return (
            solution = solution,
            functionval = fval,
            history = ithist,
            stats = stats,
            idid = idid,
            solhist = solhist,
        )
    else
        return (solution = solution, functionval = fval, 
        history = ithist, stats=stats, idid = idid)
    end
end

