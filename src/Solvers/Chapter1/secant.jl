"""
secant(f,x0; rtol=1.e-6, atol=1.e-12, maxit=10,
        armmax=10, armfix=false, pdata=nothing,
        printerr=true, keepsolhist=true, stagnationok=false)

C. T. Kelley, 2021

The secant method for scalar equations. 

Input:\n 
f: function\n 
x0: initial iterate


Keyword Arguments (kwargs):\n
rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations\n

If you use secant and your initial iterate is poor, you have made
a mistake. You will get an error message.

armmax: upper bound on stepsize reductions in linesearch

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
A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

solution = converged result
functionval = F(solution)
history = the vector of residual norms (||F(x)||) for the iteration
stats = named tuple of the history of (ifun, ijac, iarm), the number
of functions/derivatives/steplength reductions at each iteration.
For the secant method, ijac = 0.

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations
        = 1  if the line search failed

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n
secant builds solhist with a function from the Tools directory. For
systems, solhist is an N x K array where N is the length of x and K 
is the number of iteration + 1. So, for scalar equations (N=1), solhist
is a row vector. Hence the use of solhist' in the example below.


### Example for secant.jl

```jldoctest

julia> secout=secant(atan,1.0;maxit=6,atol=1.e-12,rtol=1.e-12);


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
"""
function secant(
    f,
    x0;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    solver = "secant",
    armmax = 5,
    armfix = false,
    dx = 1.e-7,
    pdata = nothing,
    printerr = true,
    keepsolhist = true,
    stagnationok = false,
)
    itc = 0
    idid = true
    errcode = 0
    iline = false
    #=
     The theory does not support convergence of the secant-Armijo iteration
     and you assume a risk when you use it. The same is true for Broyden
     and any other quasi-Newton method.
    =#
    fc = 0.0
    fc = EvalF!(f, fc, x0, pdata)
    fm = fc
    xm = copy(x0)
    xm = x0 * 1.0001
    if xm == 0
        xm = 0.0001
    end
    fm = fc
    fm = EvalF!(f, fm, xm, pdata)
    newfun0 = 1
    derivative_is_old = false
    resnorm = abs(fc)
    jfact = nothing
    stagflag = stagnationok && (armmax == 0)
    (ItRules, x, n) =
        Secantinit(x0, dx, f, solver, armmax, armfix, maxit, printerr, pdata, jfact)
    #
    # Initialize the iteration statistics
    #
    newiarm = -1
    ItData = ItStats(resnorm, 2)
    newfun = 0
    newjac = 0
    newsol = x
    xt = x
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
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
    toosoon = (resnorm <= tol)
    #
    # The main loop stops on convergence, too many iterations, or a
    # line search failure after a derivative evaluation.
    #
    while (resnorm > tol) && (itc < maxit) && (armstop || stagnationok)
        newfun = 0
        #
        # Extra function call at the start.
        #
        newjac = 0
        newfun = 0
        #
        df = (fc - fm) / (x - xm)
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
        xm = x
        x = AOUT.ax
        fm = fc
        fc = AOUT.afc
        #
        # If the line search fails and the derivative is current, 
        # stop the iteration.
        #
        armstop = AOUT.idid || derivative_is_old
        iline = ~armstop && ~stagflag
        newiarm = AOUT.aiarm
        #
        # Keep the books.
        #
        residm = resnorm
        resnorm = AOUT.resnorm
        residratio = resnorm / residm
        updateStats!(ItData, newfun, newjac, AOUT)
        #
        itc += 1
        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    (idid, errcode) = NewtonOK(resnorm, iline, tol, toosoon, itc, ItRules)
    newtonout = CloseIteration(x, fc, ItData, idid, errcode, keepsolhist, solhist)
    return newtonout
end
