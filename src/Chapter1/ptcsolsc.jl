"""
ptcsolsc(f, x, fp=difffp; rtol=1.e-6, atol=1.e-12, fp=difffp, 
        dt0=1.e-6, maxit=100, keepsolhist=true)

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

Options:\n
rtol, atol: real and absolute error tolerances\n

dt0: initial time step. The default value of 1.e-3 is a bit conservative 
and is one option you really should play with. Look at the example
where I set it to 1.0!\n

maxit: upper bound on number of nonlinear iterations. This is 
coupled to dt0. If your choice of dt0 is too small (conservative)
then you'll need many iterations to converge and will need a larger
value of maxit.

keepsolhist: if true you get the history of the iteration in the output 
tuple. This is on by default for scalar equations and off for systems.
Only turn it on if you have use for the data, which can get REALLY LARGE.

Output: A tuple (solution, functionval, history, idid, solhist) where
history is the array of absolute function values |f(x)|
of residual norms and time steps. Unless something has gone badly wrong,
dt approx |f(x_0)|/|f(x)|.

idid=true if the iteration succeeded and false if not.

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
    x,
    fp = difffp;
    rtol = 1.e-6,
    atol = 1.e-12,
    dt0 = 1.e-3,
    maxit = 100,
    keepsolhist = true,
)
    itc = 0
    idid = true
    fval = f(x)
    tol = atol + rtol * abs(fval)
    h = 1.e-7
    dt = dt0
    ithist = [abs(fval)]
    if keepsolhist
        solhist = [x]
    end
    while itc < maxit + 1 && abs(fval) > tol
        df = fpeval_newton(x, f, fval, fp, h)
        idt = 1.0 / dt
        step = -fval / (idt + df)
        x = x + step
        fvalm = fval
        fval = f(x)
        # SER 
        newhist = abs(fval)
        dt = dt * abs(fvalm) / abs(fval)
        itc = itc + 1
        if keepsolhist
            newsol = x
            append!(solhist,newsol)
        end
         append!(ithist,newhist)
    end
    resnorm=abs(fval)
    #
    if resnorm > tol
        idid = false
        PTCError(resnorm, maxit, dt0)
    end
    if keepsolhist
        return (
            solution = x,
            functionval = fval,
            history = ithist,
            idid = idid,
            solhist = solhist,
        )
    else
        return (solution = x, functionval = fval, history = ithist, idid = idid)
    end
end
