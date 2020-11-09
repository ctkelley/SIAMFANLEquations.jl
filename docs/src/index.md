# SIAMFANLEquations.jl v0.2.2

[C. T. Kelley](https://ctk.math.ncsu.edu)

[SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl) is the package of solvers and test problems
for the book

__Solving Nonlinear Equations with Iterative Methods:__
__Solvers and Examples in Julia__

__Testing github actions: Working now. The next tag is the bit test.__

This documentation is sketchy and designed to get you going, but the real deal is the [IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL)

This is version 0.2.2. 

- Registration of SIAMFANLEquations.jl
- Chapter 2 done. 

__Chapter 3__ will be Newton-Krylov solvers. That will take some time because
I need to modify my GMRES and BiCGStab solvers. I write my own Krylov
solvers so they will communicate with the nonlinear solvers the way I want
and so I can do GMRES my way (classical Gram-Schmidt twice!). The nonlinear
part won't be so bad once the linear solvers do what I want.

## Scalar Equations: Chapter 1

### Algorithms
The examples in the first chapter are scalar equations that illustrate
many of the important ideas in nonlinear solvers. 

1. infrequent reevaluation of the derivative 
2. secant equation approximation of the derivative
2. line searches
3. pseudo-transient continuation

Leaving out the kwargs, the calling sequence for getting nsolsc
to solve ``f(x) = 0`` is

```julia
nsolsc(f,x, fp=difffp)
```

Here x is the initial iterate and fp (optional) is the function
for evaluating the derivative. If you leave fp out, nsold uses
a forward difference approximation.

See the code overview or the notebook for details. Here are a couple 
of simple examples.

Solve ``atan(x) = 0`` with ``x_0 = 0`` as the initial iterate and a
finite difference approximation to the derivative. The output of
nsolsc is a tuple. The history vector contains the nonlinear residual
norm. In this example I've limited the number of iterations to 5, so
history has 6 components (including the initial residual, iteration 0).

```julia
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

Now try the same problem with the secant method. I'll need one more
iteration to meet the termination criterion.

```julia
julia> secout=secant(atan,1.0;maxit=6,atol=1.e-12,rtol=1.e-12);


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

In this example I define a function and its derivative and send that
to nsolsc. I print both the history vectors and the solution history.
```julia
julia> fs(x)=x^2-4.0; fsp(x)=2x;

julia> nsolout=nsolsc(fs,1.0,fsp; maxit=5,atol=1.e-9,rtol=1.e-9);

julia> [nsolout.solhist.-2 nsolout.history]
6×2 Array{Float64,2}:
 -1.00000e+00  3.00000e+00
  5.00000e-01  2.25000e+00
  5.00000e-02  2.02500e-01
  6.09756e-04  2.43940e-03
  9.29223e-08  3.71689e-07
  2.22045e-15  8.88178e-15
```


## Nonlinear systems with direct linear solvers: Chapter 2
The ideas from Chapter 1 remain important here. For systems the Newton step is the solution of the linear system

``F'(x) s = - F(x)``

This chapter is about solving the equation for the Newton step with Gaussian elimination. Infrequent reevaluation of ``F'``means that we also factor ``F'`` infrequently, so the impact of this idea is greater. Even better, there is typically no loss in the nonlinear iteration if we do that factorization in single precision. You an make that happen by giving nsold and ptcsold the single precision storage for the Jacobian. Half precision is also possible, but is a very, very bad idea. 

Bottom line: __single precision can cut the linear algebra cost in half with no loss in the quality of the solution or the number of nonlinear iterations it takes to get there.__


## Overview of the Codes

The core solvers (so far) are 

1. nsol.jl is is all variations of Newton's method __except__
   pseudo transient continuation. The methods are
   - Newton's method
   - The Shamanskii method, where the derivative evaluation is
     done every m iterations. ``m=1`` is Newton and ``m=\infty`` is chord.
   - I do an Armijo line search for all the methods unless the method is
     chord or you tell me not to.

2. ptcsol.jl is pseudo-transient continuation.

The solvers for scalar equations are wrappers for the core codes with
the same interface. The expectation is that the scalar codes do not need
to manage linear algebra or storage of arrays.

### Scalar Equations: Chapter 1
There are two codes for the methods in this chapter

1. nsolsc.jl is all variations of Newton's method __except__ 
   pseudo transient continuation. The methods are
   - Newton's method 
   - The Shamanskii method, where the derivative evaluation is
     done every m iterations. ``m=1`` is Newton and ``m=\infty`` is chord.
   - I do an Armijo line search for all the methods unless the method is
     chord or you tell me not to.

2. secant.jl is the scalar secant method. It is a stand-alone code.
I do not know if I'll merge it with the Broyden code or not. It's really
too simple to mess with much.

3. ptcsolsc.jl is pseudo-transient continuation. 

### Nonlinear systems with direct linear solvers: Chapter 2

This is the same story as it was for scalar equations, 'ceptin for the
linear algebra. The linear solvers for this chapter are the matrix
factorizations that live in LinearAlgebra, SuiteSparse,
or BandedMatrices. 

### Nonlinear systems with iterative linear solvers: Chapter 3

1. The plan is to make this part of nsol and ptcsol. 
2. The GMRES linear solver is 50% done. I expect to get Bi-CGSTAB in there too.


