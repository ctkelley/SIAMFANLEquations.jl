"""
    nsoli(F!, x0, FS, FPS, Jvec=dirder; rtol=1.e-6, atol=1.e-12,
               maxit=20, lmaxit=-1, lsolver="gmres", eta=.1,
               fixedeta=true, Pvec=nothing, pside="right",
               armmax=10, dx = 1.e-7, armfix=false, pdata = nothing,
               printerr = true, keepsolhist = false, stagnationok=false)
)

C. T. Kelley, 2021

Julia versions of the nonlinear solvers from my SIAM books. 
Herewith: nsoli

You must allocate storage for the function and the Krylov basis in advance
--> in the calling program <-- ie. in FS and FPS

Inputs:\n
- F!: function evaluation, the ! indicates that F! overwrites FS, your
    preallocated storage for the function.\n
    So FS=F!(FS,x) or FS=F!(FS,x,pdata) returns FS=F(x)


- x0: initial iterate\n

- FS: Preallocated storage for function. It is a vector of size N\n
  You should store it as (N) and design F! to use vectors of size (N).
  If you use (N,1) consistently instead, the solvers may work, but I make
  no guarantees.

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
    You can declare FPS as Float64 or Float32 and nsoli 
    will do the right thing. Float16 support is there, but not working well.
    
    If the Jacobian is reasonably well conditioned, you can cut the cost
    of orthogonalization and storage (for GMRES) in half with no loss. 
    There is no benefit if your linear solver is not GMRES or if 
    othogonalization and storage of the Krylov vectors is only a
    small part of the cost of the computation. So if your preconditioner
    is good and you only need a few Krylovs/Newton, reduced precision won't
    help you much.

    BiCGSTAB does not benefit from reduced precsion. 

----------------------

Keyword Arguments (kwargs):\n

rtol and atol: relative and absolute error tolerances\n

maxit: limit on nonlinear iterations\n

lmaxit: limit on linear iterations. If lmaxit > m-1, where FPS has
m columns, and you need more than m-1 linear iterations, then GMRES 
will restart. 

The default is -1 for GMRES. This means that you'll take m-1 iterations, 
where size(V) = (n,m), and get no restarts. For BiCGSTAB the default is 10.

lsolver: the linear solver, default = "gmres"\n
Your choices will be "gmres" or "bicgstab". However,
gmres is the only option for now.

eta and fixed eta: eta > 0 or there's an error

The linear solver terminates when ||F'(x)s + F(x) || <= etag || F(x) ||

where 

etag = eta if fixedeta=true

etag = Eisenstat-Walker as implemented in book if fixedeta=false

The default, which may change, is eta=.1, fixedeta=true

Pvec: Preconditioner-vector product. The rules are similar to Jvec
    So, Pv=Pvec(v,x) or Pv=Pvec(v,x,pdata) returns P(x) v where
    P(x) is the preconditioner. You must use x as an input even
    if your preconditioner does not depend on x

pside: apply preconditioner on pside, default = "right". I do not
      recommend "left". See Chapter 3 for the story on this.

armmax: upper bound on step size reductions in line search\n

dx: default = 1.e-7\n
difference increment in finite-difference derivatives
      h=dx*norm(x,Inf)+1.e-8

armfix: default = false\n
The default is a parabolic line search (ie false). Set to true and
the step size will be fixed at .5. Don't do this unless you are doing
experiments for research.\n

pdata:\n 
precomputed data for the function/Jacobian-vector/Preconditioner-vector
products.  Things will go better if you use this rather than hide the data 
in global variables within the module for your function/Jacobian

If you use pdata in any of F!, Jvec, or Pvec, you must use in in all of them.

printerr: default = true\n
I print a helpful message when the solver fails. To suppress that
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
- A named tuple (solution, functionval, history, stats, idid,
               errcode, solhist)
where

   -- solution = converged result

   -- functionval = F(solution)

   -- history = the vector of residual norms (||F(x)||) for the iteration

   -- stats = named tuple of the history of (ifun, ijac, iarm, ikfail), the 
number of functions/Jacobian-vector prods/steplength reductions/linear solver
failures at each iteration. Linear solver failures DO NOT mean that the
nonlinear solver will fail. You should look at this stat if, for example,
the line search fails. Increasing the size of FPS and/or lmaxit might
solve the problem.

I do not count the function values for a finite-difference derivative
because they count toward a Jacobian-vector product.

  -- idid=true if the iteration succeeded and false if not.

  -- errcode = 0 if if the iteration succeeded

        = -1 if the initial iterate satisfies the termination criteria

        = 10 if no convergence after maxit iterations

        = 1  if the line search failed

   -- solhist:\n
      This is the entire history of the iteration if you've set
      keepsolhist=true\n

solhist is an N x K array where N is the length of x and K is the number of
iteration + 1. So, for scalar equations, it's a row vector.

------------------------

### Example for nsoli

#### Simple 2D problem. 
You should get the same results as for nsol.jl because
GMRES will solve the equation for the step exactly in two iterations. Finite
difference Jacobians and analytic Jacobian-vector products for full precision
and finite difference Jacobian-vector products for single precision.

BiCGSTAB converges in 5 itertions and each nonlinear iteration costs
two Jacobian-vector products. Note that the storage for the Krylov
space in GMRES (jvs) is replace by a single vector (fpv) when BiCGSTAB
is the linear solver.

```jldoctest
julia> function f!(fv,x)
       fv[1]=x[1] + sin(x[2])
       fv[2]=cos(x[1]+x[2])
       return fv
       end
f! (generic function with 1 method)

julia> function JVec(v, fv, x)
       jvec=zeros(2);
       p=-sin(x[1]+x[2])
       jvec[1]=v[1]+cos(x[2])*v[2]
       jvec[2]=p*(v[1]+v[2])
       return jvec
       end
JVec (generic function with 1 method)

julia> x0=ones(2); fv=zeros(2); jv=zeros(2,2); 

julia> jv32=zeros(Float32,2,2);

julia> jvs=zeros(2,3); jvs32=zeros(Float32,2,3);

julia> nout=nsol(f!,x0,fv,jv; sham=1);

julia> kout=nsoli(f!,x0,fv,jvs,JVec; 
                  fixedeta=true, eta=.1, lmaxit=2);

julia> kout32=nsoli(f!,x0,fv,jvs32; 
                    fixedeta=true, eta=.1, lmaxit=2);

julia> [nout.history kout.history kout32.history]
5×3 Array{Float64,2}:
 1.88791e+00  1.88791e+00  1.88791e+00
 2.43119e-01  2.43120e-01  2.43119e-01
 1.19231e-02  1.19231e-02  1.19230e-02
 1.03266e-05  1.03261e-05  1.03264e-05
 1.46388e-11  1.40862e-11  1.39825e-11

julia> fpv=zeros(2);

julia> koutb=nsoli(f!,x0,fv,fpv,JVec; 
            fixedeta=true, eta=.1, lmaxit=2, lsolver="bicgstab");

julia> koutb.history
6-element Vector{Float64}:
 1.88791e+00
 2.43120e-01
 1.19231e-02
 4.87500e-04
 7.54236e-06
 3.84646e-07
```




"""
function nsoli(
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
    pside = "right",
    armmax = 10,
    dx = 1.e-7,
    armfix = false,
    pdata = nothing,
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
    stagflag = stagnationok && (armmax == 0)
    #=
    Named tuple with the iteration data. This makes communiction
    with the linear solvers and the line search easier.
    =#
    (ItRules, x, n, solhist) = Newton_Krylov_Init(
        x0,
        dx,
        F!,
        Jvec,
        Pvec,
        pside,
        lsolver,
        eta,
        fixedeta,
        armmax,
        armfix,
        maxit,
        lmaxit,
        printerr,
        pdata,
        keepsolhist,
    )
    #    keepsolhist ? (solhist = solhistinit(n, maxit, x)) : (solhist = [])
    #
    # First Evaluation of the function. Initialize the iteration stats.
    # Fix the tolerances for convergence and define the derivative FPF
    # outside of the main loop for scoping.
    #   
    FS = EvalF!(F!, FS, x, pdata)
    resnorm = norm(FS)
    tol = rtol * resnorm + atol
    FPF = []
    ItData = ItStatsK(resnorm)
    newiarm = -1
    newfun = 0
    newjac = 0
    newikfail = 0
    ke_report = false
    residratio = 1.0
    armstop = true
    etag = eta
    #
    # Preallocate a few vectors for the step, trial step, trial function
    #
    step = copy(x)
    xt = copy(x)
    FT = copy(x)
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
    toosoon = (resnorm <= tol)
    #
    # The main loop stops on convergence, too many iterations, or a
    # line search failure after a derivative evaluation.
    #
    while resnorm > tol && itc < maxit && (armstop || stagnationok)
        #   
        newfun = 0
        newjac = 0
        newikfail = 0
        #
        #
        # The GMRES solver will do the orthogonalization in lower
        # precision. I've tested Float32, but see the docstrings
        # for all the caveats. This is not the slam dunk it was
        # for Gaussian elimination on dense matrices.
        #
        step .*= 0.0
        etag = forcing(itc, residratio, etag, ItRules, tol, resnorm)
        kout = Krylov_Step!(step, x, FS, FPS, ItRules, etag)
        step .= kout.step
        #
        # For GMRES you get 1 jac-vec per iteration and there is no jac-vec
        # for the initial inner iterate of zero. For BiCGSTAB it's two 
        # jac-vecs per iteration.
        #
        newjac = kout.Lstats.lits
        (lsolver == "gmres") || (newjac *= 2)
        linok = kout.Lstats.idid
        linok || (ke_report = Krylov_Error(lmaxit, ke_report); newikfail = 1)
        #
        # Compute the trial point, evaluate F and the residual norm.     
        # The derivative is never old for Newton-Krylov
        #
        AOUT = armijosc(xt, x, FT, FS, step, resnorm, ItRules, false)
        #
        # update solution/function value
        #
        x .= AOUT.ax
        FS .= AOUT.afc
        #
        # If the line search fails 
        # stop the iteration. Print an error message unless
        # stagnationok == true
        #
        armstop = AOUT.idid
        iline = ~armstop && ~stagflag
        #
        # Keep the books.
        #
        residm = resnorm
        resnorm = AOUT.resnorm
        residratio = resnorm / residm
        updateStats!(ItData, newfun, newjac, AOUT, newikfail)
        newiarm = AOUT.aiarm
        itc += 1
        keepsolhist && (@views solhist[:, itc+1] .= x)
        #        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    #    solution = x
    #    functionval = FS
    (idid, errcode) = NewtonOK(resnorm, iline, tol, toosoon, itc, ItRules)
    newtonout = CloseIteration(x, FS, ItData, idid, errcode, keepsolhist, solhist)
    return newtonout
end
