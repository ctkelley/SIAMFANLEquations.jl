"""
nsolsc(f,x; rtol=1.e-6, atol=1.e-12, maxit=10,
        fp=difffp, solver="newton", sham=1, armmax=10, resdec=.1,
        armfix=false, printerr=true, keepsolhist=true)

Newton's method for scalar equations. Has most of the features a
code for systems of equations needs.

Input:\n 
f: function\n 
x: initial iterate\n
Options:\n
fp: derivative. If your derivative function is fp, you give me
its name. For example fp=foobar tells me that foobar is your
function for the derivative. The default is a forward difference
Jacobian that I provide.\n

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
The covergence rate has local q-order sham+1 if you only count
iteratons where you update the derivative. You need not
provide your own derivative function to use this option. sham=Inf
is chord only if chord is converging well.\n

rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations\n

sham: update Jacobian every sham iteraitons. sham=1 --> Newton

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

"""
function nsolsc(
    f,
    x;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    fp = difffp,
    solver = "newton",
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
    #
    # If you like the secant or sham=large methods, I will do a 
    # difference Jacobian anyhow if the line search kicks in. 
    # You will thank me for this.
    # Even if you don't thank me, I will do it anyhow.
    #
    # If you insist, solver=chord will ignore poor convergence and let
    # things go south with no interference. Please don't do that as
    # standard procedure.
    #
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
    #
    # Novices: If I don't define dfold outside the while loop
    # I will get a bizzare error message because of Julia's 
    # scoping rules. Don't believe me, erase that line and
    # see for yourself.
    # RTFM
    #
    dfold = 0.0
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
    while (resid > tol) && (itc < maxit)
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
                dfold = df
                derivative_is_old = false
            else
                df = dfold
                derivative_is_old = true
            end
        end
        xm = x
        fm = fc
        d = -fc / df
        AOUT = armijosc(fc, d, xm, fm, f, h, fp, armmax, 
                        armfix, derivative_is_old, printerr)
        if AOUT.idid == false
            iline = false
        end
        fc = AOUT.afc
        x = AOUT.ax
        newiarm = AOUT.aiarm
        newfun=newfun+newiarm+1
        derivative_is_old = AOUT.adfo
        d = AOUT.ad
        resid = abs(fc)
        residratio = abs(fc) / abs(fm)
        itc = itc + 1
        newhist = [abs(fc)]
        if keepsolhist
            newsol = [x]
            solhist = [solhist' newsol']'
        end
        iarm = [iarm' newiarm']'
        ifun = [ifun' newfun']'
        ijac = [ijac' newjac']'
        ithist = [ithist' newhist']'
    end
    solution = x
    fval = fc
    if abs(fval) > tol || iline == false
        idid = false
    end
    if idid == false && printerr
        println("Newton failure; maybe increase maxit and/or armmax")
        if iline == false
            println("The line search failed at least once.")
        end
        println("Current values: maxit  =  ", maxit, ", armmax = ", armmax)
        println("Give the history array a look to see what's happening.")
        println("  ")
        idid = false
    end
    stats = (ifun=ifun, ijac=ijac, iarm=iarm)
    if keepsolhist
        return (
            solution = solution,
            functionval = fval,
            history = ithist,
            stats = stats,
            iarm = iarm,
            idid = idid,
            solhist = solhist,
        )
    else
        return (solution = solution, functionval = fval, 
        history = ithist, stats=stats, idid = idid)
    end
end


"""
armijosc(fc, d, xm, fm, f, h, fp, armmax, armfix, derivative_is_old)

Line search for scalar equations. Read the notebook or print book
for the explanation. This is an internal function and I did not
design it to be hackable by the novice.
"""
function armijosc(fc, d, xm, fm, f, h, fp, armmax, armfix,
               derivative_is_old, printerr)
    idid = true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    x = xm
    lam0 = 0.0
    lamc = lambda
    lamm = lamc
    #
    #   Take the full step and, if happy, go home.
    #
    x = xm + lambda * d
    fc = f(x)
    armfail = abs(fc) > (1 - alpha * lambda) * abs(fm)
    #
    # If I have an old derivative I will not tolerate a failure in
    # the line search. 
    #
    # jflag tells me that I had to refresh the derivative
    # liarm is the counter that = iarm unless I refresh the Jacobian
    #
    jflag = false
    if armfail && derivative_is_old
        df = fpeval_newton(xm, f, fm, fp, h)
        dfold = df
        d = -fm / df
        x = xm + lambda * d
        fc = f(x)
        armfail = abs(fc) > (1 - alpha * lambda) * abs(fm)
        derivative_is_old = false
        jflag = true
    end
    liarm = 0
    iarm = 0
    #
    # At this point I've taken a full step. I'll enter the loop only if
    # that full step has failed.
    #
    ffc = abs(fc)^2
    ff0 = abs(fm)^2
    ffm = ffc
    while armfail && iarm < armmax
        #
        #   At this point fp = f'(xm) then it's time to be serious 
        #   about the line  search.
        #
        if iarm == 0 || armfix == true
            lambda = lambda * 0.5
        else
            lamm=lamc
            lamc=lambda
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm)
        end
        x = xm + lambda * d
        fc = f(x)
        ffm = ffc
        ffc = abs(fc)^2
        iarm += 1
        liarm += 1
        armfail = abs(fc) > (1 - alpha * lambda) * abs(fm)
    end
    if iarm >= armmax 
        idid = false
        if printerr
        println("Linesearch failure")
        end
    end
    return (ax = x, afc = fc, aiarm = iarm, adfo = derivative_is_old, ad = d, idid = idid)
end
