"""
ptcsol(F!, x0, FS, FPS, J! = diffjac!; rtol=1.e-6, atol=1.e-12,
             maxit=20, delta0=1.e-6, dx=1.e-7, pdata = nothing, jfact = klfact,
               printerr = true, keepsolhist = false, jknowsdt = false)

C. T. Kelley, 2021

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: some new stuff ==> ptcsol

PTC finds the steady-state solution of u' = -F(u), u(0) = u_0.
The - sign is a convention.

You must allocate storage for the function and Jacobian in advance
--> in the calling program <-- ie. in FS and FPS

Inputs:\n
- F!: function evaluation, the ! indicates that F! overwrites FS, your
    preallocated storage for the function.\n
    So, FV=F!(FV,x) or FV=F!(FV,x,pdata) returns FV=F(x)
    

- x0: initial iterate\n

- FS: Preallocated storage for function. It is a vector of size N\n
  You should store it as (N,) and design F! to use vectors of size (N,).
  If you use (N,1) consistently instead, the solvers may work, but I make
  no guarantees.

- FPS: preallocated storage for Jacobian. It is an N x N matrix\n
  If FPS is sparse, you __must__ allocate storage for the diagonal so
  I will have room to put 1/dt in there.

- J!: Jacobian evaluation, the ! indicates that J! overwrites FPS, your
    preallocated storage for the Jacobian. If you leave this out the
    default is a finite difference Jacobian.\n
    So, FP=J!(FP,FV,x) or FP=J!(FP,FV,x,pdata) returns FP=F'(x);
    (FP,FV, x) must be the argument list, even if FP does not need FV.
    One reason for this is that the finite-difference Jacobian
    does and that is the default in the solver.

    You may have a better way to add (1/dt) I to your Jacobian. If you
    want to do this yourself then your Jacobian function should be
    FP=J!(FP,FV,x,dt) or FP=J!(FP,FV,x,dt,pdata) and return
    F'(x) + (1.0/dt)*I. \n
    You will also have to set the kwarg __jknowsdt__ to true.

- Precision: Lemme tell ya 'bout precision. I designed this code for 
    full precision
    functions and linear algebra in any precision you want. You can declare
    FPS as Float64, Float32, or Float16 and ptcsol will do the right thing if
    YOU do not destroy the declaration in your J! function. I'm amazed
    that this works so easily. If the Jacobian is reasonably well
    conditioned, you can cut the cost of Jacobian factorization and
    storage in half with no loss. For large dense Jacobians and inexpensive
    functions, this is a good deal.\n
    BUT ... There is very limited support for direct sparse solvers in
    anything other than Float64. I recommend that you only use Float64
    with direct sparse solvers unless you really know what you're doing. I
    have a couple examples in the notebook, but watch out.



----------------------

Keyword Arguments (kwargs):\n

rtol and atol: relative and absolute error tolerances\n

delta0: initial pseudo time step. The default value of 1.e-3 is a bit conservative
and is one option you really should play with. Look at the example
where I set it to 1.0!\n

maxit: limit on nonlinear iterations, default=100. \n
This is coupled to delta0. If your choice of delta0 is too small (conservative)
then you'll need many iterations to converge and will need a larger
value of maxit

For PTC you'll need more iterations than for a straight-up
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

I use jfact when I call PTCUpdate to evaluate the
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
ptcsol will use backslash to compute the Newton step.

I know that this is probably not optimal in your situation, so it is 
good to pick something else, like jfact = lu.

printerr: default = true\n
I print a helpful message when the solver fails. To suppress that
message set printerr to false.

keepsolhist: default = false\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

jknowsdt: default = false\n
Set this to true if your Jacobian evaluation function retursn
F'(x) + (1/dt) I. You'll also need to follow the rules above for
the Jacobian evaluation function. I do not recommend this and if
your Jacobian is anything other than a matrix I can't promise
anything. I've tested this for matrix outputs only.

Output:\n
A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

solution = converged result
functionval = F(solution)
history = the vector of residual norms (||F(x)||) for the iteration

Unlike nsol, nsoli, or even ptcsoli, ptcsol has a fixed cost per 
iteration of one function, one Jacobian, and one Factorization. Hence
iteration statistics are not interesting and not in the output. 

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisfies the termination criteria
        = 10 if no convergence after maxit iterations

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n 
solhist is an N x K array where N is the length of x and K is the number
of iteration + 1. So, for scalar equations, it's a row vector.

### Example for ptcsol

#### The buckling beam problem. 
You'll need to use TestProblems for this to work.

```jldoctest
julia> using SIAMFANLEquations.TestProblems

julia> n=63; maxit=1000; delta = 0.01; lambda = 20.0;

julia> bdata = beaminit(n, 0.0, lambda); x = bdata.x;

julia> u0 = x .* (1.0 .- x) .* (2.0 .- x);

julia> u0 .*= exp.(-10.0 * u0);

julia> FS = copy(u0); FPS = copy(bdata.D2);

julia> pout = ptcsol( FBeam!, u0, FS, FPS, BeamJ!; 
 rtol = 1.e-10, pdata = bdata, delta0 = delta, maxit = maxit);

julia> # It takes a few iterations to get there.
       length(pout.history)
25

julia> [pout.history[1:5] pout.history[21:25]]
5×2 Array{Float64,2}:
 6.31230e+01  9.75412e-01
 7.52624e+00  8.35295e-02
 8.31545e+00  6.58797e-04
 3.15455e+01  4.12697e-08
 3.66566e+01  6.29295e-12

julia> # We get the nonnegative stedy state.
       maximum(pout.solution)
2.19086e+00
```


"""
function ptcsol(
    F!,
    x0,
    FS = [],
    FPS = [],
    J! = diffjac!;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 20,
    delta0 = 1.e-6,
    dx = 1.e-7,
    pdata = nothing,
    jfact = klfact,
    printerr = true,
    keepsolhist = false,
    jknowsdt = false,
)
    itc = 0
    idid = true
    #
    #   Initialize the iteration
    #   As with the other codes, ItRules packages all the details of
    #   the problem so it's easy to pass them around. 
    #
    (ItRules, x, n, solhist) =
        PTCinit(x0, dx, F!, J!, delta0, maxit, pdata, jfact, keepsolhist, jknowsdt)
    #
    # First Evaluation of the function. Initialize the iteration history.
    # Fix the tolerances for convergence and define the derivative FPF
    # outside of the main loop for scoping.
    #    
    FS = EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    tol = rtol * resnorm + atol
    ItData = ItStatsPTC(resnorm)
    #
    # Preallocate a vector for the step
    #
    step = copy(x)
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
    toosoon = (resnorm <= tol)
    #
    # The main loop stops on convergence or too many iterations.
    #
    delta = delta0
    while resnorm > tol && itc < maxit
        #   
        # Evaluate and factor the Jacobian; update x, F(x), and delta.  
        #
        (x, delta, FS, resnorm) = PTCUpdate(FPS, FS, x, ItRules, step, resnorm, delta)
        #
        # Keep the books
        #
        updateStats!(ItData, resnorm)
        itc += 1
        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    (idid, errcode) = PTCOK(resnorm, tol, toosoon, ItRules, printerr)
    itout = CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist)
    return (itout)
end
