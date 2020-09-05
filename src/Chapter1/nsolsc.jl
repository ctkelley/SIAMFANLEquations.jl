"""
nsolsc(f,x, fp=difffp; rtol=1.e-6, atol=1.e-12, maxit=10,
        solver="newton", sham=1, armmax=10, resdec=.1, dx=1.e-7,
        armfix=false, 
        printerr=true, keepsolhist=true, stagnationok=false)

C. T. Kelley, 2020

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
Your choices are "newton"(default), "secant", or "chord". However, 
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

dx:\n
This is the increment for forward difference, default = 1.e-7.
dx should be roughly the square root of the noise in the function.

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

stagnationok:\n
Set this to true if you want to disable the line search and either
observe divergence or stagnation. This is only useful for research
or writing a book.

Output:\n
A tuple (solution, functionval, history, stats, idid, solhist) where
history is the vector of residual norms (|f(x)|) for the iteration
and stats is a tuple of the history of (ifun, ijac, iarm), the number
of functions/derivatives/steplength reductions at each iteration.

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian evaluation. I do count them for
the secant method model.

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations
        = 1  if the line search failed

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true

Output:\n
A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

solution = converged result
functionval = F(solution)
history = the vector of residual norms (||F(x)||) for the iteration
stats = named tuple of the history of (ifun, ijac, iarm), the number
of functions/derivatives/steplength reductions at each iteration.

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian evaluation. I do count them for
the secant method model.

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations
        = 1  if the line search failed

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

# Same problem with the secant method.

```julia
julia> secout=nsolsc(atan,1.0;maxit=6,atol=1.e-12,rtol=1.e-12, solver="secant");


julia> secout.history
7-element Array{Float64,1}:
 7.85398e-01
 5.18729e-01
 5.39030e-02
 4.86125e-03
 4.28860e-06
 3.37529e-11
 2.06924e-22
```

# If you have an analytic derivative, I will use it.

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
    fp = difffp;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    solver = "newton",
    sham = 1,
    armmax = 5,
    resdec = 0.1,
    dx = 1.e-7,
    armfix = false,
    printerr = true,
    keepsolhist = true,
    stagnationok = false,
)
    itc = 0
    idid = true
    errcode=0
    iline = false
    #=
     If you like the sham=large methods, I will evaluate the derivative 
     anyhow if the line search kicks in. 
     
     The theory does not support convergence of the secant-Armijo iteration
     and you assume a risk when you use it. The same is true for Broyden
     and any other quasi-Newton method.
                    
     The chord method will ignore poor convergence and let
     things go south with no interference. Please don't do that as
     standard procedure and, if you do, don't blame me.
    =#
    fc = f(x)
    fm = fc
    xm = x
    if solver == "secant"
        xm = x * 1.0001
        if xm == 0
            xm = 0.0001
        end
        fm = f(xm)
        sham = 1
        fp = difffp
    end
    derivative_is_old = false
    resnorm = abs(fc)
    ItRules = (
        solver = solver,
        sham = sham,
        armmax = armmax,
        armfix = armfix,
        resdec = resdec,
        dx = dx,
        f = f,
        fp = fp,
    )
    #
    # Initialize the iteration statistics
    #
    newiarm = -1
    ItData = ItStats(resnorm)
    #    ItData=ItStats(history=[resnorm])
    newfun = 0
    newjac = 0
    newsol = x
    xt = x
    if keepsolhist
        solhist = [x]
    end
    #
    # Fix the tolerances for convergence and define the derivative df
    # outside of the main loop for scoping.
    #
    tol = rtol * resnorm + atol
    residratio = 1
    df = 0.0
    armstop = true
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
    toosoon = false
    resnorm > tol || (toosoon = true)
    #
    # The main loop stops on convergence, too many iterations, or a
    # line search failure after a derivative evaluation.
    #
    while (resnorm > tol) && (itc < maxit) && (armstop || stagnationok)
        newfun = 0
        newjac = 0
        #
        # Evaluate the derivativce if (1) you are using the secant method, 
        # (2) you are using the chord method and it's the intial iterate, or
        # (3) it's Newton and you are on the right part of the Shamaskii loop,
        # or the line search failed with a stale deriviative, or the residual
        # reduction ratio is too large. This logic is a bit tedious, so I
        # put it in a function. See src/Tools/test_evaljac.jl
        #
        #        evaljacit = (itc % sham == 0 || newiarm > 0 || residratio > resdec)
        #        chordinit = (solver == "chord") && itc == 0
        #        evaljac = (evaljacit && solver == "newton") || chordinit ||
        #            solver == "secant"
        evaljac = test_evaljac(ItRules, itc, newiarm, residratio)
        # 
        # We've evaluated a derivative if the solver is Newton or we just
        # initialized the chord method. For secant it costs an extra function.
        #
        if evaljac
            df = PrepareDerivative(ItRules, x, xm, fc, fm)
            newfun += solver == "secant"
            newjac += ~(solver == "secant")
        end
        derivative_is_old = (newjac == 0) && (solver == "newton")
        #
        # Compute the Newton direction and call the line search.
        #
        xm = x
        fm = fc
        ft = fc
        d = -fc / df
        AOUT = armijosc(xt, x, ft, fc, d, resnorm, ItRules, derivative_is_old)
        #
        # update solution/function value
        #
        x = AOUT.ax
        fc = AOUT.afc
        #
        # If the line search fails and the derivative is current, 
        # stop the iteration.
        #
        armstop = AOUT.idid || derivative_is_old
        newiarm = AOUT.aiarm
        iline = ~armstop
        #
        # Keep the books.
        #
        residm = resnorm
        resnorm = AOUT.resnorm
        residratio = resnorm / residm
        updateStats!(ItData, newfun, newjac, AOUT)
        #
        itc += 1
        if keepsolhist
            newsol = x
            append!(solhist, newsol)
        end
    end
    solution = x
    fval = fc
    resfail = (resnorm > tol)
    idid = ~(resfail || iline || toosoon)
    if ~idid 
        errcode = NewtonError(resfail, iline, resnorm, toosoon, tol,
                    itc, maxit, armmax,printerr)
    end
    stats = (ifun = ItData.ifun, ijac = ItData.ijac, iarm = ItData.iarm)
    if keepsolhist
        return (
            solution = solution,
            functionval = fval,
            history = ItData.history,
            stats = stats,
            idid = idid,
            errcode = errcode,
            solhist = solhist,
        )
    else
        return (
            solution = solution,
            functionval = fval,
            history = ItData.history,
            stats = stats,
            idid = idid,
            errcode = errcode
        )
    end
end