# SIAMFANLEquations.jl

[C. T. Kelley](https://ctk.math.ncsu.edu)

[SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl) is the package of solvers and test problems
for the book

__Solving Nonlinear Equations with Iterative Methods:__
__Solvers and Examples in Julia__

## Scalar Equations: Chapter 1

### Algorithms
The examples in the first chapter are scalar equations that illustrate
many of the important ideas in nonlinear solvers. 

1. infrequent reevaluation of the derivative 
2. secant equation approximation of the derivative
2. line searches
3. pseudo-transient continuation

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

2. ptcsc.jl is pseudo-transient continuation. 

## Nonlinear systems with direct linear solvers: Chapter 2

## Nonlinear systems with iterative linear solvers: Chapter 3
