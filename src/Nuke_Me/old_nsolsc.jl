"""
nsolsc(x,f; rtol=1.e-6, atol=1.e-12, maxit=10,
        fp=difffp, solver="newton", sham=1, armmax=10, armfix=false,
        keepsolhist=true)

Newton's method for scalar equations. Has most of the features a
code for systems of equations needs.

Input:\n 
x: initial iterate\n
f: function\n 
Options:\n
fp: derivative. If your derivative function is fp, you give me
its name. For example fp=foobar tells me that foobar is your
function for the derivative. The default is a forward difference
Jacobian that I provide.\n

Your best bet is to put f and fp in the same file. See the atan
example.

solver:\n
Your choices are "newton"(default) or "secant". However, you have sham
at your disposal only if you chose newton.\n

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

I only turn Shamanskii on if the residuals are decreasing
rapidly, at least a factor of 10, and the line search is quiescent.\n  

rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations\n

sham: update Jacobian every sham iteraitons. sham=1 --> Newton

armmax: upper bound on stepsize reductions in linesearch

armfix:\n
The default is a parabolic line search (ie false). Set to true and
the stepsize will be fixed at .5. Don't do this unless you are doing
experiments for research.

keepsolhist:\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

Output:\n
A tuple (solution, functionval, history, idid, solhist) where
history is a tuple (iteration counter, f(x), iarm) with the history
of the entire iteration. iarm is the counter for steplength reductions.

idid=true if the iteration succeeded and false if not.

solhist=entire history of the iteration if keepsolhist=true
"""
function nsolsc(
    x,
    f;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    fp = difffp,
    solver = "newton",
    sham = 1,
    armmax = 5,
    armrule="constant",
    keepsolhist=true
)
    itc = 0
    idid=true
    iline=true
    h = 1.e-7
    #
    # If you like the secant or chord methods, I will do a difference Jacobian
    # anyhow if the line search kicks in. You will thank me for this.
    # Even if you don't thank me, I will do it anyhow.
    #
    if solver == "secant"
        xm = x * 1.0001
        if xm == 0
            xm = .0001
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
    iarm = 0
    ithist = [itc fc iarm]
    if keepsolhist
       solhist=[x]
    end
    tol = rtol * resid + atol
    residratio = 1
    while (resid > tol) && (itc < maxit)
        if solver == "secant"
            df = (fc - fm) / (x - xm)
        else
            if itc % sham == 0 || iarm > 0 || residratio > .1
                df = fpeval_newton(x, f, fc, fp, h)
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
        iarm = -1
        AOUT = armijo(fc, d, xm, fm, f, h, fp, armmax, derivative_is_old)
        if AOUT.idid == false
            iline=false
        end
        fc = AOUT.afc
        x = AOUT.ax
        iarm = AOUT.aiarm
        derivative_is_old = AOUT.adfo
        d = AOUT.ad
        resid = abs(fc)
        residratio = abs(fc) / abs(fm)
        itc = itc + 1
        newhist = [itc fc iarm]
        if keepsolhist
           newsol=[x]
           solhist=[solhist' newsol']'
        end
        ithist = [ithist' newhist']'
    end
    solution = x
    fval = fc
if abs(fval) > tol || iline==false
    idid=false
end
if idid==false
    println("Newton failure; maybe increase maxit and/or armmax")
    if iline==false
       println("The line search failed at least once.")
    end
    println("Current values: maxit  =  ",maxit,", armmax = ",armmax)
    println("Give the history array a look to see what's happening.")
    println("  ")
    idid=false
end
if keepsolhist
    return (solution=solution, functionval=fval, history=ithist, idid=idid,
        solhist=solhist)
    else
    return (solution=solution, functionval=fval, history=ithist, idid=idid)
end
end

#function fpeval_newton(x, f, fc, fp, h)
#    if fp == difffp
#        df = difffp(x, f, fc, h)
#    else
#        df = fp(x)
#    end
#    return df
#end

#function difffp(x, f, fc, h)
#    df = (f(x + h) - fc) / h
#    return df
#end

function armijo(fc, d, xm, fm, f, h, fp, armmax, derivative_is_old)
    idid=true
    alpha = 1.e-4
    iarm = -1
    lambda = 1.0
    x = xm
    lam0=0.0
    lamc=lambda
    lamm=lamc
    armfail=abs(fc) > (1 - alpha * lambda) * abs(fm)
    
#
#   jflag tells me that I had to refresh the derivative
#   liarm is the counter that = iarm unless I refresh the Jacobian
#
    jflag=false
    liarm=-1
    while armfail && iarm < armmax
        #
        # If I have an old derivative I will not tolerate a failure in
        # the line search. 
        #
        if iarm == 0 && derivative_is_old
            df = fpeval_newton(xm, f, fm, fp, h)
            dfold = df
            d = -fm / df
            lambda = 1.0
            derivative_is_old = false
            jflag=true
            liarm=-1
        end
        #
        #   If fp = f'(xm) then it's time to be serious about the line
        #   search.
        #
        x = xm + lambda * d
        fc = f(x)
        iarm += 1
        liarm += 1
        armfail=abs(fc) > (1 - alpha * lambda) * abs(fm)
        lambda = lambda * .5
    end
    if iarm >= armmax
       idid=false
       println("Linesearch failure")
    end
    return (ax = x, afc = fc, aiarm = iarm, adfo = derivative_is_old, ad = d,
            idid=idid)
end

