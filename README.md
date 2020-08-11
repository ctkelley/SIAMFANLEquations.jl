| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] |


# SIAMFANLEquations version 0.2.1

This is the package with the solvers and test problems for 

# Solving Nonlinear Equations with Iterative Methods: <br> Solvers and Examples in Julia

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This will be a sequel to my book 

(Kel03) C. T. Kelley, [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003.

Hence the notebook and this package all have SIAMFANL in their names.

The new book with have a different algorithm mix and the solvers and examples will be in Juila. The project will have three parts.

   1. A print book (details coming sooner or later)
   2. [An IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL) (open source)
   3. This package (MIT License)<br>
      I will register this package no earlier that the release of v0.2.1 and no later than the release of v0.2.2 (when Chapter 2 of the notebook and the solvers/test problems for that chapter are closer to done than they are now). 
   
## Readme Contents:
- [Installation](#Installation)
- [Meaning of Version Numbers](#Meaning-of-version-numbers)
- [__Please__ No Pull Requests](#Pull-Requests)
- [Core References and Documentation](#Core-References-and-Documentation)
- [Algorithms and Solvers](#Algorithms-and-Solvers)
- [About the test problems](#Test-Problems)
- [Funding](#Funding) 


## Installation: Use Julia 1.5 and up with this thing!!!

This package has been tested on Julia 1.4 and 1.5. __It no longer works on 1.0!__ Before mid-September __I will make changes that only 1.5 supports.__
You have been warned. 

Type this in the REPL to install

```
] add https://github.com/ctkelley/SIAMFANLEquations.jl
```

then, as usual
```
using SIAMFANLequations
```
enables you to use the codes. You'll need
```
using SIAMFANLEquations.TestProblems
```
to run the test problems.



## Meaning of version numbers

If __log(version_number) < 0__ there's trouble!

This is the dev version of v0.2.1. __Chapter 2: Systems of equations with direct linear solvers.__ 

I have added new keyword arguments to nsolsc.jl and make the error messages more informative. None of this should break the notebooks for Chapter 1. I will be updating the notebooks soon. 

__I have changed the user interface and the calling sequence of the solvers
since v0.1.1 and will continue to do things like that until v1.0.0 goes out the door.__

I have released version v0.1.2. The codes can now duplicate the examples (at least the ones that will remain in the new book)  from Chapter 1 of (Kel03), make all the new examples, and I have finished Chaper 1 of the notebook and print book. The plan is, for x > 1

-- v0.x.1 goes live when the codes can duplicate the examples I'll keep from Chapter x of (Kel03) and make the new examples. Version v0.x.2 goes out when the codes and notebook for Chapter x are finished. 0.x.y for y > 2 are serious bug fixes and/or changes in the calling sequences/interfaces/rules that I have to do to make things consistent with future chapters.

-- I will formally register the package no earlier than version 0.2.1 and no later than version 0.2.2. 

-- v1.0.0 goes out when the print book is published. This means that after v1.0.0 the interface to the codes will always be consistent with the book. My readers get my __solemn word__ on that.

## Pull Requests

__Please, please__, do not send me PRs. If you find a bug in the codes or an error in the documentation/notebook, please 

- tell me the old fashioned way with email to tim_kelley@ncsu.edu
- or open an issue.

This is a book project and I need to put all changes in by hand so I'll have muscle memory about what's going on. If you object to an algorithmic choice, you'll have to be content to know that I have thought about the algorithm mix pretty carefully and understand this field fairly well.

## Core References and Documentation

The best documentation for this package will be the [notebook](https://github.com/ctkelley/NotebookSIAMFANL) and the print book. They will have detailed algorithmic descriptions, examples for you to play with, and guidance on tweaking the algorithmic paramenters to solve your problems. The notebook will be built in parallel with the print book and the content will be __roughly__ the same. The differences will be to accommodate the two formats. For example, docstrings need some work after the map from notebook to print and notebook has to make sense as an interactive resource.

I've also used [documenter.jl](https://github.com/JuliaDocs/Documenter.jl) with this package. Click the badge
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/stable)
to get the documentation from the latest release. The documenter files have the headers for the solvers and some of the test problems. I continue to work on the docs and they will get better, but will never be as good as the notebook.

This book will not cover theory in detail (ie no proofs). My two books on nonlinear equations

(Kel95) C. T. Kelley, [***Iterative Methods for Linear and Nonlinear Equations***](https://my.siam.org/Store/Product/viewproduct/?ProductId=862) , Frontiers in Applied Mathematics 16,  SIAM 1995

and

(Kel03) C. T. Kelley, [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003

describe the Newton and Broyden algorithms. Kel95 has the theory. This project is a sequal to Kel03. Kel03 is Matlab-centric
and will remain in print.

A recent Acta Numerica paper has everything

(Kel18) C. T. Kelley, ***Numerical Methods for Nonlinear Equations***, Acta Numerica 27 (2018), pp 207--287. https://doi.org/10.1017/S0962492917000113

The references I use for theory of pseudo-transient continuation and Anderson acceleration are

(KK98) C. T. Kelley and D. E. Keyes, ***Convergence Analysis of Pseudo-Transient Continuation***, SIAM Journal on Numerical Analysis 35 (1998), pp 508-523. https://doi.org/10.1137/S0036142996304796

(TK15) A. Toth and C. T. Kelley, ***Convergence Analysis for Anderson Acceleration***, SIAM Journal on Numerical Analysis 53, (2015), pp 805-819. https://doi.org/10.1137/130919398

## Algorithms and Solvers

The solvers are designed to be stand-alone codes. The reason for this is the education mission of the project. I want the codes to be as easy to understand as possible. I have deliberately sacrificed a lot of abstraction and some performance in this effort. The reward for the reader is that the algorithmic parameters are completely exposed so  you can play with them. At the end I hope to write a wrapper for all this that hides the parameters, but the stand-alone, keyword-infested codes are what you need if you want to really understand how these methods work. My students became experts in this field by fiddling with the Matlab version of these solvers.

The linear solvers are tuned to communicate well with nonlinear solvers. My old Matlab codes are a good illustration of this idea. My [new Mablab codes](https://ctk.math.ncsu.edu/knl.html) were designed in response to the need to do this better than I had been. In particular, the linear solver and the matrix-vector/preconditioner-vector product function need information on the nonlinear iteration and any precomputed data. While I could use global variables (and did in Kel95) and put these things in a module to simplify the interface, I won't do that anymore if I can avoid it. Global varaibles break parallelism and I like to avoid them. 

The algorithms, listed by book chapter will be

   - Chapter 1: Newton-Armijo and Pseudo-transient continuation for scalar equations: nsolsc.jl and ptcsc.jl
        - Codes: __Done!__, Notebook: __Done!__
   - Chapter 2: Newton-Armijo and Pseudo-transient continuation for systems with direct linear solvers: nsold.jl and ptcd.jl
        - Codes: __In progress. The examples with dense Jacobians work.__
   - Chapter 3: Newton-Armijo and Pseudo-transient continuation for systems with iterative linear solvers: nsoli.jl and ptci.jl
       - Linear solver(s): klgmres.jl and maybe klbicgstab.jl
   - Chapter 4: Anderson acceleration: aasol.jl
   - Chapter 5: Broyden's method: brsol.jl
   
   
## Test Problems

You'll need the TestProblems submodule to run the notebook. To get it type __using SIAMFANLEquations.TestProblems__ in the repl or the notebook.

There are two kinds of test problems. The ones you care about are the ones that I use in the print book and notebook to demonstrate the algorithms. The "inside baseball" problems are the ones I __only__ use for CI. They only appear in the /test directory. If you don't know or care what CI is, be happy.
   
## Funding

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

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ctkelley.github.io/SIAMFANLEquations.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ctkelley.github.io/SIAMFANLEquations.jl/stable

[travis-img]: https://travis-ci.com/ctkelley/SIAMFANLEquations.jl.svg?branch=master
[travis-url]: https://travis-ci.com/ctkelley/SIAMFANLEquations.jl

[codecov-img]: https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl
