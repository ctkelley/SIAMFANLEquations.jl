"""
ptcsold(F!, x0, FS, FPS, J! = diffjac!; rtol=1.e-6, atol=1.e-12,
               maxit=20, dt0=1.e-6, dx=1.e-7, pdata = nothing, jfact = klfact,
               printerr = true, keepsolhist = false)

C. T. Kelley, 2020

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: ptcsold

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
    FPS as Float64, Float32, or Float16 and ptcsold will do the right thing if
    YOU do not destroy the declaration in your J! function. I'm amazed
    that this works so easily. If the Jacobian is reasonably well 
    conditioned, I can see no reason to do linear algebra in 
    double precision. 

----------------------

Keyword Arguments (kwargs):\n

rtol and atol: relative and absolute error tolerances\n

dt0: initial time step. The default value of 1.e-3 is a bit conservative
and is one option you really should play with. Look at the example
where I set it to 1.0!\n

maxit: limit on nonlinear iterations, default=100. \n
This is coupled to dt0. If your choice of dt0 is too small (conservative)
then you'll need many iterations to converge and will need a larger
value of maxit

For PTC you'll need more iterations than for a stright-up
nonlinear solve. This is part of the price for finding the 
stable solution.

dx: default = 1.e-7\n
difference increment in finite-difference derivatives
      h=dx*norm(x)+1.e-6

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
ptcsold will use backslash to compute the Newton step.

I know that this is probably not optimal in your situation, so it is 
good to pick something else, like jfact = lu.

Please do not mess with the line that calls PrepareJac!. 
        FPF = PrepareJac!(FPS, FS, x, ItRules,dt)
PrePareJac! builds F'(u) + (1/dt)I, factors it, and sends ptcsold that
factorization.

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

"""
function ptcsold(
    F!,
    x0,
    FS=[],
    FPS=[],
    J! = diffjac!;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 20,
    dt0 = 1.e-6,
    dx = 1.e-7,
    pdata = nothing,
    jfact = klfact,
    printerr = true,
    keepsolhist = false,
)
    itc = 0
    idid = true
    iline = false
    #
    #    First evaluation of the function. 
    #
    x = zeros(size(x0))
    n = length(x0)
    x .= x0
    if keepsolhist
        solhist = zeros(n, maxit + 1)
        @views solhist[:, 1] .= x
    end
    EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    ithist = [resnorm]
    residm = resnorm
    ItRules = (dx = dx, f = F!, fp = J!, pdata = pdata, fact = jfact)

    #
    # Initialize the iteration statistics
    #   
    newfun = 0
    newjac = 0
    #
    # Fix the tolerances for convergence and define the derivative FPF
    # outside of the main loop for scoping.
    #    
    tol = rtol * resnorm + atol
    FPF = []
    #
    # Preallocate a few vectors for the step, trial step, trial function
    #
    step = zeros(size(x))
    xt = zeros(size(x))
    FT = zeros(size(x))
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
    toosoon = false
    resnorm > tol || (toosoon = true)
    #
    # The main loop stops on convergence, too many iterations, or a
    # line search failure after a derivative evaluation.
    #
    dt = dt0
    while resnorm > tol && itc < maxit
        #   
        # Evaluate and factor the Jacobian.   
        #
        newfun = 0
        newjac = 0
        #
        FPF = PrepareJac!(FPS, FS, x, ItRules, dt)
        #
        step .= -(FPF \ FS)
        #
        # update solution/function value
        #
        x .= x + step
        EvalF!(F!, FS, x, pdata)
        resnorm = norm(FS)
        append!(ithist, resnorm)
        #
        # Update dt
        #
        dt *= (residm / resnorm)
        #
        # Keep the books
        #
        residm = resnorm
        itc += 1
        if keepsolhist
            @views solhist[:, itc+1] .= x
        end
    end
    solution = x
    functionval = FS
    resfail = (resnorm > tol)
    idid = ~(resfail || toosoon)
    errcode = 0
    if ~idid
        errcode = PTCError(resnorm, maxit, dt0, toosoon, tol, printerr)
    end
    if keepsolhist
        sizehist = itc + 1
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            idid = idid,
            errcode = errcode,
            solhist = solhist[:, 1:sizehist],
        )
    else
        return (
            solution = x,
            functionval = FS,
            history = ithist,
            idid = idid,
            errcode = errcode,
        )
    end
end
