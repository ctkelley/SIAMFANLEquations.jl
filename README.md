[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/dev)
[![Build Status](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl.svg?branch=master)](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl)
[![codecov](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl)
# SIAMFANLEquations

This is the package with the solvers and test problems for 

# Solving Nonlinear Equations with Iterative Methods: <br> Solvers and Examples in Julia

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This will be a sequal to my book 

(Kel03) [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003.

Hence the notebook and this package all have SIAMFANL in their names.

The new book with have a different algorithm mix and the solvers and examples will be in Juila. The project will have three parts.

   1. A print book (details coming sooner or later)
   2. [An IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL) (open source)
   3. This package (MIT License)<br>
      I will register this package when Chapter 2 of the notebook and the solvers/test problems it needs are done.
   
## Readme Contents:

- [Core References and Documentation](#Core-References-and-Documentation)
- [Algorithms and Solvers](#Algorithms-and-Solvers)
- [Support](#Support)


## Core References and Documentation

The best documentation for this pacakge will be the [notebook](https://github.com/ctkelley/NotebookSIAMFANL). It will have detailed algorithmic descriptions, examples for you to play with, and guidance on tweaking the algorithmic paramenters to solve your problems. The notebook will be built in parallel project with the print book and the content will be roughly the same.

I've also used [documenter.jl[(https://github.com/JuliaDocs/Documenter.jl) with this package. Click the badge
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/dev)
to get the headers for the solvers and some of the test problems. I continue to work on the docs and they will get better, but will never be as good as the notebook.

This book will not cover theory. My two books on nonlinear equations

(Kel95) [***Iterative Methods for Linear and Nonlinear Equations***](https://my.siam.org/Store/Product/viewproduct/?ProductId=862) , Frontiers in Applied Mathematics 16,  SIAM 1995

and

(Kel03) [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003

describe the Newton and Broyden algoirthms. CTK95 has the theory. This project is a sequal to CTK03. CTK03 is Matlab-centric
and will remain in print.

A recent Acta Numerica paper has everything

(Kel18) C. T. Kelley, ***Numerical Methods for Nonlinear Equations***, Acta Numerica 27 (2018), pp 207--287. https://doi.org/10.1017/S0962492917000113

The references I use for theory of pseudo-transient continuation and Anderson acceleration are

(KK98) C. T. Kelley and D. E. Keyes, ***Convergence Analysis of Pseudo-Transient Continuation***, SIAM Journal on Numerical Analysis 35 (1998), pp 508-523. https://doi.org/10.1137/S0036142996304796

(TK15) A. Toth and C. T. Kelley, ***Convergence Analysis for Anderson Acceleration***, SIAM Journal on Numerical Analysis 53, (2015), pp 805-819. https://doi.org/10.1137/130919398

## Algorithms and Solvers

The solvers are designed to be stand-alone codes. The reason for this is the education mission of the project. I want the codes to be as easy to understand as possible. I have deliberately sacrificed a lot of abstraction and some performance in this effort. The reward for the reader is that the algorithmic parameters are completely exposed so  you can play with them. At the end I plan to write a wrapper for all this that hides the parameters, but the stand-alone, keyword-infested codes are what you need if you want to really understand how these methods work. My students became experts in this field by fiddling with the Matlab version of these solvers.

The linear solvers are tuned to communicate well with nonlinear solvers. My old Matlab codes are a good illustration of this idea. My [new Mablab codes](https://ctk.math.ncsu.edu/knl.html) were designed in response to the need to do this better than I had been. In particular, the linear solver and the matrix-vector/preconditioner-vector product function need information on the nonlinear iteration and any precomputed data. While I could use global variables and put these things in a module to simplify the interface, I won't. Global varaibles break parallelism and I like to avoid them. I have had to use globals a couple times and am trying to get rid of them.

The algorithms, listed by book chapter will be

   - Chapter 1: Newton-Armijo and Pseudo-transient continuation for scalar equations: nsolsc.jl and ptcsc.jl
   - Chapter 2: Newton-Armijo and Pseudo-transient continuation for systems with direct linear solvers: nsold.jl and ptcd.jl
   - Chapter 3: Newton-Armijo and Pseudo-transient continuation for systems with iterative linear solvers: nsoli.jl and ptci.jl
   - Chapter 4: Anderson acceleration: anderson.jl
   - Chapter 5: Broyden's method: brsol.jl
   - Chapter 6: Linear solver(s): klgmres.jl and maybe klbicgstab.jl
   
## Support

This project was partially supported by
1. Army Research Office grant W911NF-16-1-0504 and
2. National Science Foundation Grants
   1. OAC-1740309
   2. DMS-1745654
   3. DMS-1906446
   
Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author and
do not necessarily reflect the views of the National
Science Foundation
or the Army Research Office.
