"""
nsolsc(f,x0, fp=difffp; rtol=1.e-6, atol=1.e-12, maxit=10,
        solver="newton", sham=1, armmax=10, resdec=.1, dx=1.e-7,
        armfix=false, 
        printerr=true, keepsolhist=true, stagnationok=false)

C. T. Kelley, 2020

Newton's method for scalar equations. Has most of the features a
code for systems of equations needs.

Input:\n 
f: function\n 
x0: initial iterate\n
fp: derivative. If your derivative function is fp, you give me
its name. For example fp=foobar tells me that foobar is your
function for the derivative. The default is a forward difference
Jacobian that I provide.\n


Keyword Arguments (kwargs):\n
rtol, atol: real and absolute error tolerances\n

maxit: upper bound on number of nonlinear iterations\n

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
The convergence rate has local q-order sham+1 if you only count
iterations where you update the derivative. You need not
provide your own derivative function to use this option. sham=Inf
is chord only if chord is converging well.\n

armmax: upper bound on stepsize reductions in linesearch

resdec: target value for residual reduction. \n
The default value is .1. In the old MATLAB codes it was .5.
I only turn Shamanskii on if the residuals are decreasing
rapidly, at least a factor of resdec, and the line search is quiescent.
If you want to eliminate resdec from the method ( you don't ) then set
resdec = 1.0 and you will never hear from it again.  

dx:\n
This is the increment for forward difference, default = 1.e-7.
dx should be roughly the square root of the noise in the function.

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

stagnationok:\n
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
because they count toward a Jacobian evaluation. I do count them for
the secant method model.

idid=true if the iteration succeeded and false if not.

errcode = 0 if if the iteration succeeded
        = -1 if the initial iterate satisifies the termination criteria
        = 10 if no convergence after maxit iterations
        = 1  if the line search failed

solhist:\n
This is the entire history of the iteration if you've set
keepsolhist=true\n
nsolsc builds solhist with a function from the Tools directory. For
systems, solhist is an N x K array where N is the length of x and K 
is the number of iteration + 1. So, for scalar equations (N=1), solhist
is a row vector. Hence the use of solhist' in the example below.


# Examples
```jldoctest
julia> nsolout=nsolsc(atan,1.0;maxit=5,atol=1.e-12,rtol=1.e-12);

julia> nsolout.history
6-element Array{Float64,1}:
 7.85398e-01
 5.18669e-01
 1.16332e-01
 1.06102e-03
 7.96200e-10
 2.79173e-24
```

# Same problem with the secant method.

```julia
julia> secout=nsolsc(atan,1.0;maxit=6,atol=1.e-12,rtol=1.e-12, solver="secant");


julia> secout.history
7-element Array{Float64,1}:
 7.85398e-01
 5.18729e-01
 5.39030e-02
 4.86125e-03
 4.28860e-06
 3.37529e-11
 2.06924e-22
```

# If you have an analytic derivative, I will use it.

```jldoctest
julia> fs(x)=x^2-4.0; fsp(x)=2x;

julia> nsolout=nsolsc(fs,1.0,fsp; maxit=5,atol=1.e-9,rtol=1.e-9);

julia> [nsolout.solhist'.-2 nsolout.history]
6×2 Array{Float64,2}:
 -1.00000e+00  3.00000e+00
  5.00000e-01  2.25000e+00
  5.00000e-02  2.02500e-01
  6.09756e-04  2.43940e-03
  9.29223e-08  3.71689e-07
  2.22045e-15  8.88178e-15

```

"""
function nsolsc(
    f,
    x0,
    fp = difffp;
    rtol = 1.e-6,
    atol = 1.e-12,
    maxit = 10,
    solver = "newton",
    sham = 1,
    armmax = 5,
    resdec = 0.1,
    dx = 1.e-7,
    armfix = false,
    printerr = true,
    keepsolhist = true,
    stagnationok = false,
)
    itc = 0
    idid = true
    errcode=0
    iline = false
    #=
     If you like the sham=large methods, I will evaluate the derivative 
     anyhow if the line search kicks in. 
     
     The theory does not support convergence of the secant-Armijo iteration
     and you assume a risk when you use it. The same is true for Broyden
     and any other quasi-Newton method.
                    
     The chord method will ignore poor convergence and let
     things go south with no interference. Please don't do that as
     standard procedure and, if you do, don't blame me.
    =#
    fc = f(x0)
    fm = fc
    xm = copy(x0)
    if solver == "secant"
        xm = x0 * 1.0001
        if xm == 0
            xm = 0.0001
        end
        fm = f(xm)
        sham = 1
        fp = difffp
    end
   fp0=copy(x0)
   fpp0=copy(x0)
   zdata=[]
   newtonout=nsol(f, x0, fp0, fpp0, fp;
         rtol=rtol,atol=atol,maxit=maxit,solver=solver, sham=sham,
         armmax=armmax, resdec=resdec, dx=dx, armfix=armfix,
       pdata=zdata,printerr=printerr,keepsolhist=keepsolhist, 
       stagnationok=stagnationok)
    return newtonout
end
