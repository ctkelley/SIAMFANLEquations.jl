"""
    nsold(F!, x0, FS, FPS, J!=diffjac!; rtol=1.e-6, atol=1.e-12,
               maxit=20, solver="newton", sham=1, armmax=10, resdec=.1,
               dx = 1.e-7, armfix=false, 
               pdata = nothing, jfact = lu!,
               printerr = true, keepsolhist = false, stagnationok=false)
)

C. T. Kelley, 2020

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: nsold


Inputs:\n
- F!: function evaluation, the ! indicates that F! overwrites FS, your
    preallocated storage for the function.\n

- x: initial iterate\n

- FS: Preallcoated storage for function. It is an N x 1 column vector\n

- FPS: preallcoated storage for Jacobian. It is an N x N matrix\n

- J!: Jacobian evaluation, the ! indicates that J! overwrites FPS, your
    preallocated storage for the Jacobian. If you leave this out the
    default is a finite difference Jacobian.\n

----------------------

Keyword Arguments (kwargs):\n

- rtol and atol: relative and absolute error tolerances\n

- maxit: limit on nonlinear iterations\n

solver:\n
Your choices are "newton"(default) or "chord". However,
you have sham at your disposal only if you chose newton. "chord"
will keep using the initial derivative until the iterate converges,
uses the iteration budget, or the line search fails. It is not the
same as sham=Inf, which is smarter.\n

sham:\n
This is the Shamanskii method. If sham=1, you have Newton.
The iteration updates the derivative every sham iterations.
The convergence rate has local q-order sham+1 if you only count
iterations where you update the derivative. You need not
provide your own derivative function to use this option. sham=Inf
is chord only if chord is converging well.\n

armmax: upper bound on stepsize reductions in linesearch\n

resdec:\n 
target value for residual reduction.
The default value is .1. In the old MATLAB codes it was .5.
I only turn Shamanskii on if the residuals are decreasing
rapidly, at least a factor of resdec, and the line search is quiescent.
If you want to eliminate resdec from the method ( you don't ) then set
resdec = 1.0 and you will never hear from it again.

dx:\n
difference increment in finite-difference derivatives
      h=dx*norm(x)+1.e-6

armfix:\n
The default is a parabolic line search (ie false). Set to true and
the stepsize will be fixed at .5. Don't do this unless you are doing
experiments for research.\n

pdata:\n 
precomputed data for the function/Jacobian. 
Things will go better if you use this rather than hide the data 
in global variables within the module for your function/Jacobian

jfact:\n
If you have a dense Jacobian I call PrepareJac! to evaluate the
Jacobian (using your J!) and factor it. The default is to use
lu! to compute an LU factorization and share storage with the
Jacobian. You may change LU to something else by, for example,
setting jfact = cholseky! if your Jacobian is spd. 

Please do not mess with the line that calls PrepareJac!. 
        FPF=PrepareJac!(FS, FPS, x, F!, J!, dx, pdata; fact = jfact)
FPF is not the same as FPS (the storage you allocate for the Jacobian)
for a reason. FPF and FPS do not have the same type, even though they
share storage. So, FPS=PrepareJac!(FS, FPS, ...) will break things.

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

------------------------

# Using nsold.jl

Here are the rules as of August 11, 2020.

F! is the nonlinear residual. 
J! is the Jacobian evaluation.

I like to put all my function/Jacobian/initialization stuff in a
Module and only export the things I actually use.

A) You allocate storage for the function and Jacobian in advance 
   --> in the calling program <-- NOT in FS and FPS

FV=F!(FV,x) or FV=F!(FV,x,pdata) returns FV=F(x)

FP=J!(FV,FP,x) or FP=J!(FV,FP,x,pdata) returns FP=F'(x); 
    (FV,FP, x) must be the argument list, even if FP does not need FV.
    One reason for this is that the finite-difference Jacobian
    does and that is the default in the solver.

In the future J! will also be a matrix-vector product and FPS will
be the PREALLOCATED (!!) storage for the GMRES(m) Krylov vectors.

Lemme tell ya 'bout precision. I designed this code for full precision
functions and linear algebra in any precision you want. You can declare
FPS as Float64, Float32, or Float16 and nsold will do the right thing if 
YOU do not destroy the declaration in your J! function. I'm amazed 
that this works so easily. 

If the Jacobian is reasonably well conditioned, I can see no reason
to do linear algebra in double precision

Don't try to evaluate function and Jacobian all at once because 
that will cost you a extra function evaluation every time the line
search kicks in.

B) Any precomputed data for functions, Jacobians, matrix-vector products
   or preallocted storage may live in global variables within a module 
   containing F! and J!.  Don't do that if you can avoid it. 
   Use pdata instead.

# Examples
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
help?> heqinit
search: heqinit

  heqinit(x0::Array{T,1}, c, TJ=Float64) where T<:Real

  Initialize H-equation precomputed data.

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
    jfact = lu!,
    printerr = true,
    keepsolhist = false,
    stagnationok = false,
)

    itc = 0
    idid = true
    iline = true
    #=
    First evaluation of the function. I evaluate the derivative when
    Shamanskii tells me to, at the first iteration (duh!), and when
    the rate of residual reduction is below the target value of resdec.
    =#
    x = zeros(size(x0))
    n = length(x0)
    x .= x0
    if keepsolhist
        solhist = zeros(n, maxit + 1)
        @views solhist[:, 1] .= x
    end
    EvalF!(FS, x, F!, pdata)
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
    while resnorm > tol && itc < maxit && (armstop || stagnationok)
        #   
        # Evaluate and factor the Jacobian.   
        #
        newfun = 0
        newjac = 0
        #
        # Evaluate the derivativce if (1) you are using the chord method 
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
            FPF = PrepareJac!(FS, FPS, x, ItRules)
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
        # stop the iteration.
        #
        armstop = AOUT.idid || derivative_is_old
        iline = ~armstop
        #
        # Keep the books.
        #
        residm = resnorm
        resnorm = AOUT.resnorm
        residratio = resnorm / residm
        updateStats!(ItData, newfun, newjac, AOUT)
        newiarm = AOUT.aiarm
        itc += 1
        if keepsolhist
            @views solhist[:, itc+1] .= x
        end
    end
    solution = x
    functionval = FS
    resfail = (resnorm > tol)
    idid = ~(resfail || iline)
    if ~idid && printerr
        NewtonError(resfail, iline, resnorm, itc, maxit, armmax)
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
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ItData.history,
            stats = stats,
            idid = idid,
        )
    end
end
