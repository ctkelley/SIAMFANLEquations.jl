"""
function ptcsoli(
    F!,
    x0,
    FS,
    FPS,
    Jvec = dirder;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 20,
    lmaxit = -1,
    lsolver = "gmres",
    eta = 0.1,
    fixedeta = true,
    Pvec = nothing,
    PvecKnowsdelta = false, 
    pside = "right",
    delta0 = 1.e-6,
    dx = 1.e-7,
    pdata = nothing,
    printerr = true,
    keepsolhist = false,
)

C. T. Kelley, 2021

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: some new stuff ==> ptcsoli

PTC finds the steady-state solution of u' = -F(u), u(0) = u_0.
The - sign is a convention.

You must allocate storage for the function and Krylov basis in advance
--> in the calling program <-- ie. in FS and FPS

Inputs:\n
- F!: function evaluation, the ! indicates that F! overwrites FS, your
    preallocated storage for the function.\n
    So, FV=F!(FV,x) or FV=F!(FV,x,pdata) returns FV=F(x)

- x0: initial iterate\n

- FS: Preallocated storage for function. It is an N x 1 column vector.\n
You may dimension it as (n,) or (n,1). (n,) is best, but the
solvers can deal with it either way.

- FPS: preallocated storage for the Krylov basis. It is an N x m matrix where
       you plan to take at most m-1 GMRES iterations before a restart. \n

- Jvec: Jacobian vector product, If you leave this out the
    default is a finite difference directional derivative.\n
    So, FP=Jvec(v,FS,x) or FP=Jvec(v,FS,x,pdata) returns FP=F'(x) v. \n
    (v, FS, x) or (v, FS, x, pdata) must be the argument list,
    even if FP does not need FS.
    One reason for this is that the finite-difference derivative
    does and that is the default in the solver.

- Precision: Lemme tell ya 'bout precision. I designed this code for 
    full precision functions and linear algebra in any precision you want. 
    You can declare FPS as Float64 or Float32 and ptcsoli 
    will do the right thing. Float16 support is there, but not working well.

    If the Jacobian is reasonably well conditioned, you can cut the cost
    of orthogonalization and storage (for GMRES) in half with no loss.
    There is no benefit if your linear solver is not GMRES or if
    othogonalization and storage of the Krylov vectors is only a
    small part of the cost of the computation. So if your preconditioner
    is good and you only need a few Krylovs/Newton, reduced precision won't
    help you much.

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
stable solution. \n

lmaxit: limit on linear iterations. If lmaxit > m-1, where FPS has
m columns, and you need more than m-1 linear iterations, then GMRES
will restart.

The default is -1. For GMRES this means that you'll take m-1 iterations, where
size(V) = (n,m), and get no restarts. For BiCGSTAB you'll then get the default
of 10 iterations.

lsolver: the linear solver, default = "gmres"\n
Your choices will be "gmres" or "bicgstab". However,
gmres is the only option for now. \n

eta and fixed eta: eta > 0 or there's an error.

The linear solver terminates when ||F'(x)s + F(x) || <= etag || F(x) ||

where

etag = eta if fixedeta=true

etag = Eisenstat-Walker as implemented in book if fixedeta=false

The default, which may change, is eta=.1, fixedeta=true \n

Pvec: Preconditioner-vector product. The rules are similar to Jvec
    So, Pv=Pvec(v,x) or Pv=Pvec(v,x,pdata) returns P(x) v where
    P(x) is the preconditioner. You must use x as an input even
    if your preconditioner does not depend on x.\n 

PvecKnowsdelta: If you want your preconditioner-vector product to depend on 
    the pseudo-timestep delta, put an array deltaval in your precomputed
    data. Initialize it as
    deltaval = zeros(1,)
    and let ptcsoli know about it by setting the kwarg
    PvecKnowsdelta = true
    ptcsoli will update the value in deltaval with every change
    to delta with pdata.deltaval[1]=delta
    so your preconditioner-vector product can get to it.\n
    

pside: apply preconditioner on pside, default = "right". I do not
      recommend "left". The problem with "left" for ptcsoli is
      that it can fail to satisfy the inexact Newton condition for 
      the unpreconditioned equation, especially early in the iteration
      and lead to an incorrect result (unstable solution or wrong 
      branch of steady state).
      See Chapter 3 for the story on this. \n

dx: default = 1.e-7\n
difference increment in finite-difference derivatives
      h=dx*norm(x)+1.e-8 \n

pdata:\n 
precomputed data for the function/Jacobian-vector/Preconditioner-vector
products.  Things will go better if you use this rather than hide the data
in global variables within the module for your function/Jacobian

If you use pdata in any of F!, Jvec, or Pvec, you must use in in all of them.
precomputed data for the function/Jacobian. 
Things will go better if you use this rather than hide the data 
in global variables within the module for your function/Jacobian. \n

printerr: default = true\n
I print a helpful message when the solver fails. To suppress that
message set printerr to false. \n

keepsolhist: default = false\n
Set this to true to get the history of the iteration in the output
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.\n

Output:\n
A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

solution = converged result
functionval = F(solution)
history = the vector of residual norms (||F(x)||) for the iteration
stats = named tuple of the history of (ifun, ijac, ikfail), the number
of functions/jacobian-vector prodcuts/linear solver filures at each iteration.

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian-vector product.

Linear solver failures need not cause the nonlinear iteration to fail. 
You get a warning and that is all. \n

idid=true if the iteration succeeded and false if not. \n

errcode = 0 if if the iteration succeeded \n
        = -1 if the initial iterate satisfies the termination criteria
        = 10 if no convergence after maxit iterations \n

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n 
solhist is an N x K array where N is the length of x and K is the number
of iteration + 1. So, for scalar equations, it's a row vector.

## Example from the docstrings for ptcsol

### The buckling beam problem. 
You'll need to use TestProblems for
this to work. The preconditioner is a solver for the high order term.

```jldoctest
julia> using SIAMFANLEquations.TestProblems

julia> function PreCondBeam(v, x, bdata)
          J = bdata.D2
          ptv = J\v
       end
PreCondBeam (generic function with 1 method)

julia> n=63; maxit=1000; delta0 = 0.01; lambda = 20.0;

julia> bdata = beaminit(n, 0.0, lambda);

julia> x = bdata.x; u0 = x .* (1.0 .- x) .* (2.0 .- x); u0 .*= exp.(-10.0 * u0);


julia> FS = copy(u0); FPJV=zeros(n,20);

julia> pout = ptcsoli( FBeam!, u0, FS, FPJV; delta0 = delta0, pdata = bdata,
       eta = 1.e-2, rtol = 1.e-10, maxit = maxit, Pvec = PreCondBeam);

julia> # It takes a few iterations to get there.
       length(pout.history)
25

julia> [pout.history[1:5] pout.history[21:25]]
5×2 Array{Float64,2}:
 6.31230e+01  1.79578e+00
 7.45926e+00  2.65964e-01
 8.73598e+00  6.58278e-03
 2.91936e+01  8.35069e-06
 3.47969e+01  5.11594e-09

julia> # We get the nonnegative stedy state.
       norm(pout.solution,Inf)
2.19086e+00

n=63; maxit=1000; delta0 = 0.01; lambda = 20.0;

julia> # Use BiCGSTAB for the linear solver

julia> FS = copy(u0); FPJV=zeros(n,);

julia> pout = ptcsoli( FBeam!, u0, FS, FPJV; delta0 = delta0, pdata = bdata,
       eta = 1.e-2, rtol = 1.e-10, maxit = maxit, 
       Pvec = PreCondBeam, lsolver="bicgstab");

julia> # Same number of iterations as GMRES, but each one costs double 

julia> # the Jacobian-vector products and much less storage

julia> length(pout.history)
25

julia> [pout.history[1:5] pout.history[21:25]]
5×2 Matrix{Float64}:
 6.31230e+01  1.68032e+00
 7.47081e+00  2.35073e-01
 8.62095e+00  5.18262e-03
 2.96495e+01  3.23715e-06
 3.51504e+01  3.33107e-10

```

"""
function ptcsoli(
    F!,
    x0,
    FS,
    FPS,
    Jvec = dirder;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 20,
    lmaxit = -1,
    lsolver = "gmres",
    eta = 0.1,
    fixedeta = true,
    Pvec = nothing,
    PvecKnowsdelta = false,
    pside = "right",
    delta0 = 1.e-6,
    dx = 1.e-7,
    pdata = nothing,
    printerr = true,
    keepsolhist = false,
)
    itc = 0
    idid = true
    #
    #   Initialize the iteration
    #   As with the other codes, ItRules packages all the details of
    #   the problem so it's easy to pass them around. 
    #
    (ItRules, x, n) = PTCKrylovinit(
        x0,
        dx,
        F!,
        Jvec,
        delta0,
        Pvec,
        PvecKnowsdelta,
        pside,
        lsolver,
        eta,
        fixedeta,
        lmaxit,
        maxit,
        printerr,
        pdata,
    )
    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    #
    # First Evaluation of the function. Initialize the iteration history.
    # Fix the tolerances for convergence and define the derivative FPF
    # outside of the main loop for scoping.
    #    
    FS = EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    ItData = ItStatsPTCK(resnorm)
    tol = rtol * resnorm + atol
    etag = eta
    ke_report = false
    residratio = 1.0
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
        residm = resnorm
        newjac = 0
        newikfail = 0
        #   
        # Comppute the Jacobian-vector product; update x, F(x), and delta.  
        #
        etag = forcing(itc, residratio, etag, ItRules, tol, resnorm)
        (x, delta, FS, resnorm, Lstats) =
            PTCUpdatei(FPS, FS, x, ItRules, step, resnorm, delta, etag)
        resdiratio = resnorm / residm
        newjac = Lstats.lits
        linok = Lstats.idid
        linok || (ke_report = Krylov_Error(lmaxit, ke_report); newikfail = 1)
        #
        # Keep the books
        #
        updateStats!(ItData, resnorm, newjac, newikfail)
        itc += 1
        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    (idid, errcode) = PTCOK(resnorm, tol, toosoon, ItRules, printerr)
    itout = CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist)
    return itout
end
