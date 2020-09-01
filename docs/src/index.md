# SIAMFANLEquations.jl v0.2.1

[C. T. Kelley](https://ctk.math.ncsu.edu)

[SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl) is the package of solvers and test problems
for the book

__Solving Nonlinear Equations with Iterative Methods:__
__Solvers and Examples in Julia__

This documentation is sketchy and designed to get you going, but the real deal is the [IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL)

This is version 0.2.1. 

** This thing is under constant revision. I think the user interfaces to
nsolsc, ptcsolsc, and nsold are stable, but you never know.***

The scalar solvers and the first chapter of the notebook are done as of
v0.1.2.

Chapter 2 is under construction and I'll tag this when the solvers are 
done. I'll tag v0.2.2 when the notebook is ready.

nsold.jl, Newton with direct linear solvers, is done. I am finishing 
the test problems now.

The notebooks for Chapter 2 are nowhere close to done. The to-do list inlcudes

0. Finishing the test problems and the solvers. (75% done)
1. Getting the print book part of Chapter 2 looking the way I want. (75% done)
2. Making the formatting of Chapter 1 consistent with Chapter 2. (25% done)
3. Mapping the print book part of Chapter 2 to the notebook. (0% done)
4. Completing the notebook part of Chapter 2. (10% done)
5. Mapping the notebook part of Chapter 2 to the printbook. (0% done)

If all goes well, I should post a draft of everything by late September

Once item 0 is done I will tag v0.2.1.

## Scalar Equations: Chapter 1

### Algorithms
The examples in the first chapter are scalar equations that illustrate
many of the important ideas in nonlinear solvers. 

1. infrequent reevaluation of the derivative 
2. secant equation approximation of the derivative
2. line searches
3. pseudo-transient continuation

See the code overview or the notebook for details. Here are a couple 
of simple examples.

Solve ``atan(x) = 0`` with ``x_0 = 0`` as the initial iterate and a
finite difference approximation to the derivative.

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

## Nonlinear systems with direct linear solvers: Chapter 2
The ideas from Chapter 1 remain important here. For systems the Newton step is the solution of the linear system

``F'(x) s = - F(x)``

This chapter is about solving the equation for the Newton step with Gaussian elimination. Infrequent reevaluation of ``F'``means that we also factor ``F'`` infrequenly, so the impact of this idea is greater. Even better, there is typically no loss in the nonlinear iteration if we do that factorization in single precision. You an make that happen by giving nsold and ptcsold the single precision storage for the Jacobian. Half precision is also possible, but is a very, very bad idea. 

Bottom line: __single precision can cut the linear algebra cost in half with no loss in the quality of the solution or the number of nonlinear iterations it takes to get there.

## Nonlinear systems with iterative linear solvers: Chapter 3

## Overview of the Codes

### Scalar Equations: Chapter 1
There are two codes for the methods in this chapter

1. nsolsc.jl is all variations of Newton's method __except__ 
   pseudo transient continuation. The methods are
   - Newton's method 
   - The Shamanskii method, where the derivative evaluation is
     done every m iterations. ``m=1`` is Newton and ``m=\infty`` is chord.
   - The secant method
   - You have the option to do an Armijo line search for all the methods.

2. ptcsolsc.jl is pseudo-transient continuation. 

### Nonlinear systems with direct linear solvers: Chapter 2

### Nonlinear systems with iterative linear solvers: Chapter 3
