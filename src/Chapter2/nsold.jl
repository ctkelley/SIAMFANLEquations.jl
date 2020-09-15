"""
    nsold(F!, x0, FS, FPS, J!=diffjac!; rtol=1.e-6, atol=1.e-12,
               maxit=20, solver="newton", sham=1, armmax=10, resdec=.1,
               dx = 1.e-7, armfix=false, 
               pdata = nothing, jfact = klfact,
               printerr = true, keepsolhist = false, stagnationok=false)
)

C. T. Kelley, 2020

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: nsold

You must allocate storage for the function and Jacobian in advance
--> in the calling program <-- ie. in FS and FPS

Inputs:\n
- F!: function evaluation, the ! indicates that F! overwrites FS, your
    preallocated storage for the function.\n

    FV=F!(FV,x) or FV=F!(FV,x,pdata) returns FV=F(x)

- x0: initial iterate\n

- FS: Preallcoated storage for function. It is an N x 1 column vector\n

- FPS: preallcoated storage for Jacobian. It is an N x N matrix\n

- J!: Jacobian evaluation, the ! indicates that J! overwrites FPS, your
    preallocated storage for the Jacobian. If you leave this out the
    default is a finite difference Jacobian.\n

    FP=J!(FP,FV,x) or FP=J!(FP,FV,x,pdata) returns FP=F'(x);
    (FP,FV, x) must be the argument list, even if FP does not need FV.
    One reason for this is that the finite-difference Jacobian
    does and that is the default in the solver.

    Lemme tell ya 'bout precision. I designed this code for full precision
    functions and linear algebra in any precision you want. You can declare
    FPS as Float64, Float32, or Float16 and nsold will do the right thing if
    YOU do not destroy the declaration in your J! function. I'm amazed
    that this works so easily. If the Jacobian is reasonably well 
    conditioned, I can see no reason to do linear algebra in 
    double precision. 

----------------------

Keyword Arguments (kwargs):\n

- rtol and atol: relative and absolute error tolerances\n

- maxit: limit on nonlinear iterations\n

solver: default = "newton"\n
Your choices are "newton" or "chord". However,
you have sham at your disposal only if you chose newton. "chord"
will keep using the initial derivative until the iterate converges,
uses the iteration budget, or the line search fails. It is not the
same as sham=Inf, which is smarter.\n

sham: default = 1 (ie Newton)\n
This is the Shamanskii method. If sham=1, you have Newton.
The iteration updates the derivative every sham iterations.
The convergence rate has local q-order sham+1 if you only count
iterations where you update the derivative. You need not
provide your own derivative function to use this option. sham=Inf
is chord only if chord is converging well.\n

armmax: upper bound on stepsize reductions in linesearch\n

resdec: default = .1\n 
This is the target value for residual reduction.
The default value is .1. In the old MATLAB codes it was .5.
I only turn Shamanskii on if the residuals are decreasing
rapidly, at least a factor of resdec, and the line search is quiescent.
If you want to eliminate resdec from the method ( you don't ) then set
resdec = 1.0 and you will never hear from it again.

dx: default = 1.e-7\n
difference increment in finite-difference derivatives
      h=dx*norm(x)+1.e-6

armfix: default = false\n
The default is a parabolic line search (ie false). Set to true and
the stepsize will be fixed at .5. Don't do this unless you are doing
experiments for research.\n

pdata:\n 
precomputed data for the function/Jacobian. 
Things will go better if you use this rather than hide the data 
in global variables within the module for your function/Jacobian

jfact: default = klfact (tries to figure out best choice) \n
If your Jacobian has any special structure, please set jfact
to the correct choice for a factorization.

I use jfact when I call PrepareJac! to evaluate the
Jacobian (using your J!) and factor it. The default is to use
klfact (an internal function) to do something reasonable.
For general matrices, klfact picks lu! to compute an LU factorization
and share storage with the Jacobian.  You may change LU to something else by,
for example, setting jfact = cholseky! if your Jacobian is spd.

klfact knows about banded matrices and picks qr. You should,
however RTFM, allocate the extra two upper bands, and use jfact=qr!
to override klfact.

If you give me something that klfact does not know how to dispatch on,
then nothing happens. I just return the original Jacobian matrix and 
nsold will use backslash to compute the Newton step.

I know that this is probably not optimal in your situation, so it is 
good to pick something else, like jfact = lu.

Please do not mess with the line that calls PrepareJac!. 
        FPF = PrepareJac!(FPS, FS, x, ItRules)
FPF is not the same as FPS (the storage you allocate for the Jacobian)
for a reason. FPF and FPS do not have the same type, even though they
share storage. So, FPS=PrepareJac!(FPS, FS, ...) will break things.

printerr: default = true\n
I print a helpful message when the solver fails. To supress that
message set printerr to false.

keepsolhist: default = false\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

stagnationok: default = false\n
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

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian evaluation. 

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations
        = 1  if the line search failed

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true

------------------------

# Examples
## World's easiest problem example.

```jldoctest
 julia> function f!(fv,x)
       fv[1]=x[1] + sin(x[2])
       fv[2]=cos(x[1]+x[2])
       end
f (generic function with 1 method)

julia> x=ones(2,); fv=zeros(2,); jv=zeros(2,2);
julia> nout=nsold(f!,x,fv,jv);
julia> nout.history
5-element Array{Float64,1}:
 1.88791e+00
 2.43119e-01
 1.19231e-02
 1.03266e-05
 1.46416e-11

julia> nout.solution
2-element Array{Float64,1}:
 -7.39085e-01
  2.30988e+00

```

## H-equation example

```jldoctest
julia> n=16; x0=ones(n,); FV=ones(n,); JV=ones(n,n);
julia> hdata=heqinit(x0, .5);
julia> hout=nsold(heqf!,x0,FV,JV;pdata=hdata);
julia> hout.history
3-element Array{Float64,1}:
 6.17376e-01
 3.17810e-03
 6.22034e-08
```
"""
function nsold(
    F!,
    x0,
    FS,
    FPS,
    J! = diffjac!;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 20,
    solver = "newton",
    sham = 1,
    armmax = 10,
    resdec = 0.1,
    dx = 1.e-7,
    armfix = false,
    pdata = nothing,
    jfact = klfact,
    printerr = true,
    keepsolhist = false,
    stagnationok = false,
)

    itc = 0
    idid = true
    iline = false
#
#   If I'm letting the iteration stagnate and turning off the
#   linesearch, then the line search cannot fail.
#
    stagflag = stagnationok && (armmax==0)
    #=
    First evaluation of the function. I evaluate the derivative when
    Shamanskii tells me to, at the first iteration (duh!), and when
    the rate of residual reduction is below the target value of resdec.
    =#
    x = zeros(size(x0))
    n = length(x0)
    x .= x0
    ~keepsolhist || (solhist=solhistinit(n, maxit, x))
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    ItRules = (
        solver = solver,
        sham = sham,
        armmax = armmax,
        armfix = armfix,
        resdec = resdec,
        dx = dx,
        f = F!,
        fp = J!,
        pdata = pdata,
        fact = jfact,
    )

    #
    # Initialize the iteration statistics
    #   
    newiarm = -1
    #    ItData=ItStats(history=[resnorm])
    ItData = ItStats(resnorm)
    newfun = 0
    newjac = 0
    #
    # Fix the tolerances for convergence and define the derivative FPF
    # outside of the main loop for scoping.
    #    
    tol = rtol * resnorm + atol
    derivative_is_old = false
    residratio = 1.0
    FPF = []
    armstop = true
    #
    # Preallocate a few vectors for the step, trial step, trial function
    #
    step = zeros(size(x))
    xt = zeros(size(x))
    FT = zeros(size(x))
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
#    toosoon = false
#    resnorm > tol || (toosoon = true)
    toosoon = (resnorm <= tol)
    #
    # The main loop stops on convergence, too many iterations, or a
    # line search failure after a derivative evaluation.
    #
    while resnorm > tol && itc < maxit && (armstop || stagnationok)
        #   
        # Evaluate and factor the Jacobian.   
        #
        newfun = 0
        newjac = 0
        #
        # Evaluate the derivative if (1) you are using the chord method 
        # and it's the intial iterate, or
        # (2) it's Newton and you are on the right part of the Shamaskii loop,
        # or the line search failed with a stale deriviative, or the residual
        # reduction ratio is too large. This leads to a tedious barrage
        # of conditionals that I have parked in a function.
        #
        #        evaljacit = (itc % sham == 0 || newiarm > 0 || residratio > resdec)
        #        chordinit = (solver == "chord") && itc == 0
        #        evaljac = test_evaljac(itc, solver, sham, newiarm, residratio, resdec)
        evaljac = test_evaljac(ItRules, itc, newiarm, residratio)
        if evaljac
            FPF = PrepareJac!(FPS, FS, x, ItRules)
            newjac += 1
        end
        derivative_is_old = (newjac == 0) && (solver == "newton")
        #        derivative_is_old = ~evaljacit && (solver == "newton")
        step .= -(FPF \ FS)
        #
        # Compute the trial point, evaluate F and the residual norm.     
        #
        AOUT = armijosc(xt, x, FT, FS, step, resnorm, ItRules, derivative_is_old)
        #
        # update solution/function value
        #
        x .= AOUT.ax
        FS .= AOUT.afc
        #
        # If the line search fails and the derivative is current,
        # stop the iteration. Print an error message unless
        # stagnationok == true and armmax=0
        #
        armstop = AOUT.idid || derivative_is_old
        iline = ~armstop && ~stagflag
        #
        # Keep the books.
        #
        residm = resnorm
        resnorm = AOUT.resnorm
        residratio = resnorm / residm
        updateStats!(ItData, newfun, newjac, AOUT)
        newiarm = AOUT.aiarm
        itc += 1
        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    solution = x
    functionval = FS
    resfail = (resnorm > tol)
    idid = ~(resfail || iline || toosoon)
    errcode=0
    if ~idid 
        errcode= NewtonError(resfail, iline, resnorm, toosoon, tol,
                        itc, maxit, armmax, printerr)
    end
    stats = (ifun = ItData.ifun, ijac = ItData.ijac, iarm = ItData.iarm)
    if keepsolhist
        sizehist = itc + 1
        return (
            solution = x,
            functionval = FS,
            history = ItData.history,
            stats = stats,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ItData.history,
            stats = stats,
            idid = idid,
            errcode = errcode
        )
    end
end
