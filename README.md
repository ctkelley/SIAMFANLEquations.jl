[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/dev)
[![Build Status](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl.svg?branch=master)](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl)
[![codecov](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl)
# SIAMFANLEquations

This is the package with the solvers and test problems for 

# Solving Nonlinear Equations with Iterative Methods: <br> Solvers and Examples in Julia

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This will be a sequal to my book 

(Kel03) [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003

The new book with have a different algorithm mix and the solvers and examples will bein Juila. The project will have three parts.

   1. A print book (details coming sooner or later)
   2. [An IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL) (open source)
   3. This package (MIT License)
   
## Readme Contents:


[Algorithms and Solvers](#Algorithms-and-Solvers)

## Algorithms

The algorithms, listed by book chapter will be

   - Chapter 1: Newton-Armijo and Pseudo-transient continuation for scalar equations: nsolsc.jl and ptcsc.jl
   - Chapter 2: Newton-Armijo and Pseudo-transient continuation for systems with direct linear solvers: nsold.jl and ptcd.jl
   - Chapter 3: Newton-Armijo and Pseudo-transient continuation for systems with iterative linear solvers: nsoli.jl and ptci.jl
