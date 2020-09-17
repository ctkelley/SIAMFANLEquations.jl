"""
ptcsolsc(f, x0, fp=difffp; rtol=1.e-6, atol=1.e-12, maxit=100,
        dt0=1.e-6, dx=1.e-7, printerr = true, keepsolhist=true)

C. T. Kelley, 2020

Scalar pseudo-transient continuation solver. PTC is designed to find
stable steady state solutions of 

dx/dt = - f(x)

It is ABSOLUTELY NOT a general purpose nonlinear solver.

Input:\n
f: function\n
x: initial iterate/data\n
fp: derivative. If your derivative function is fp, you give me
its name. For example fp=foobar tells me that foobar is your
function for the derivative. The default is a forward difference
Jacobian that I provide.\n

Keyword Arguments:\n
rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations. This is 
coupled to dt0. If your choice of dt0 is too small (conservative)
then you'll need many iterations to converge and will need a larger
value of maxit.

dt0: initial time step. The default value of 1.e-3 is a bit conservative 
and is one option you really should play with. Look at the example
where I set it to 1.0!\n

dx: default = 1.e-7\n
difference increment in finite-difference derivatives
      h=dx*norm(x)+1.e-6

printerr: default = true\n
I print a helpful message when the solver fails. To supress that
message set printerr to false.

keepsolhist: if true you get the history of the iteration in the output 
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

Output: A tuple (solution, functionval, history, idid, errcode, solhist) where
history is the array of absolute function values |f(x)|
of residual norms and time steps. Unless something has gone badly wrong,
dt approx |f(x_0)|/|f(x)|.

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations

solhist=entire history of the iteration if keepsolhist=true

If the iteration fails it's time to play with the tolerances, dt0, and maxit.
You are certain to fail if there is no stable solution to the equation.

# Examples
```jldoctest
julia> ptcout=ptcsolsc(sptest,.2;dt0=2.0,rtol=1.e-3,atol=1.e-3);

julia> [ptcout.solhist ptcout.history]
7×2 Array{Float64,2}:
 2.00000e-01  9.20000e-02
 9.66666e-01  4.19962e-01
 8.75086e-01  2.32577e-01
 7.99114e-01  1.10743e-01
 7.44225e-01  4.00926e-02
 7.15163e-01  8.19395e-03
 7.07568e-01  4.61523e-04
```

"""
function ptcsolsc(
    f,
    x0,
    fp = difffp;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 100,
    dt0 = 1.e-3,
    dx = 1.e-7,
    printerr = true,
    keepsolhist = true
)
    itc = 0
    idid = true
    fval=1.0
    (ItRules, x, n)=PTCinit(x0, dx, f, fp, nothing, nothing)
    keepsolhist ?  (solhist=solhistinit(n, maxit, x)) : (solhist=[])
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
#    fval = f(x)
    fval=EvalF!(f, fval, x)
    resnorm=abs(fval)
    ithist = [resnorm]
    tol = atol + rtol * resnorm
    #
    # If the initial iterate satisfies the termination criteria, tell me.
    #
    toosoon = (resnorm <= tol)
    df=0.0
    step=0.0
    dt = dt0
    while itc < maxit && resnorm > tol
#
# Take a PTC step: update the point and dt
#
(x, dt, fval, resnorm) = PTCUpdate(df, fval, x, ItRules, step, resnorm, dt)
#
# Keep the books
#
        append!(ithist, resnorm)
        itc += 1
        ~keepsolhist || (@views solhist[:, itc+1] .= x)
    end
    #
    errcode = 0
    resfail = (resnorm > tol)
    idid = ~(resfail || toosoon)
    errcode = 0
    ~(resfail || toosoon) || 
     (errcode = PTCError(resnorm, maxit, dt0, toosoon, tol, printerr))
    itout=PTCClose(x, fval, ithist, idid, errcode, keepsolhist, solhist)
    return(itout)
end
