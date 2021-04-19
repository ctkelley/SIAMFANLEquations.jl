| **Documentation**                                                               | **Build Status**                                                                                | **DOI**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------- |
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][build-status-img]][build-status-url] [![][codecov-img]][codecov-url] | [![DOI](https://zenodo.org/badge/256312455.svg)](https://zenodo.org/badge/latestdoi/256312455) |


# SIAMFANLEquations version 0.3.1
[changelog](#Changes)


This is the package with the solvers and test problems for 

# Solving Nonlinear Equations with Iterative Methods: <br> Solvers and Examples in Julia

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This will be a sequel to my book 

(Kel03) C. T. Kelley, [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM 2003.

Hence the notebook and this package all have SIAMFANL in their names.

The new book with have a different algorithm mix and the solvers and examples will be in Juila. The project will have three parts.

   1. A print book: __Under contract with SIAM for manuscript delivery in 2021 and publication in 2022__.
   2. [An IJulia notebook](https://github.com/ctkelley/NotebookSIAMFANL/releases/tag/v0.2.3) (open source, MIT License, Creative Commons License)
      - If you're using the Chapter 3 stuff in the notebook you must the [IJulia notebook Master Branch](https://github.com/ctkelley/NotebookSIAMFANL). It has the new stuff and is not likely to break things in the older version of the package __today__.
   3. This package (MIT License)<br>
   
## Readme Contents:
- [Mission](#Package-Mission)
- [Installation](#Installation)
- [Meaning of Version Numbers](#Meaning-of-version-numbers)
- [__Please__ No Pull Requests](#Pull-Requests)
- [Core References and Documentation](#Core-References-and-Documentation)
- [Algorithms and Solvers](#Algorithms-and-Solvers)
- [About the test problems](#Test-Problems)
- [How to cite this stuff](#Citations)
- [What's new in this version since 0.2.3](#Changes)
- [Funding](#Funding) 

## Package Mission

This package is designed and built to support a book project. So the solvers and examples reinforce the algorithmic discussion in the book. General purpose packages have a different mission.

## Installation: Use Julia 1.5 and up with this thing!!!

This package has been tested on Julia 1.5. __It no longer works on 1.0!__ It may still work on 1.4, but I make no promises.

Type this 

```
] add SIAMFANLEquations
```

or this

```
import Pkg; Pkg.add("SIAMFANLEquations")
```
in the REPL to install the package.

Then, as usual
```
using SIAMFANLequations
```
enables you to use the codes. You'll need
```
using SIAMFANLEquations.TestProblems
```
to run the test problems. Then there the examples you get with
```
using SIAMFANLEquations.Examples
```
for the unit tests, the examples in the book, and the notebook.


## Meaning of version numbers

If __log(version_number) < 0__ there's trouble!

This is version v0.3.1. I will tag this once PTC-Krylo
      
- nsol.jl, ptcsol.jl (Newton and pseudo-transient continuation) codes are stable. The scalar codes nsolsc.jl and ptcsol.jl are also stable.
     
The plan is, for x > 2.

- v0.x.0 goes live when the codes can duplicate the examples I'll keep from Chapter x of (Kel03) and make the new examples. 

- Version v0.x.1 means the codes are finished and I have solid drafts of the print book and notebook parts of the chapter. 

- Version v0.x.2 goes out when the codes and notebook for Chapter x are finished.

- v0.x.3 is reserved for finalizing the print book <--> notebook mappings, cleaing up the docs, and fixing inconsitencies. I post the package announcements for v0.x.3.

 - 0.x.y for y > 3 are serious bug fixes and/or changes in the calling sequences/interfaces/rules that I have to do to make things consistent with future chapters. 

v1.0.0 goes out __when the print book is published__. This means that after v1.0.0 the interface to the codes will always be consistent with the book. My readers get my __solemn word__ on that.

## Pull Requests

__I like bug reports; I need bug reports__, but ...

__Please, please__, do not send me PRs. If you find 
   1.  a bug (programming or performance) in the codes,
   2. confusion, lack of clarity, or __errors in the installation instructions__,
       1. I would __really like__ some Windows users to try this stuff, especially the notebooks.
   3. something I could do in the user interface to help you do your work ...
       1. that won't break other stuff, 
       2. make the code or __user interface__ opaque to a novice,
       3. or eat up lots of time,
   4. a factual error in the documentation/notebook, or 
   5. an error/inconsistency in the docstrings, please ...
  
 Do your choice of ... 

- tell me the old fashioned way with email to tim_kelley@ncsu.edu 
- or open an issue.

This is a book project and I need to put all changes in by hand so I'll have muscle memory about what's going on.

I have limited bandwidth, __so please do not send me email or open issues about__ ...

   1. Typos in the notebook or the docstrings. This project is far from the final proofreading stage and I want to fix those things in peace. There are many of them and I do not need 100s of emails/issues about that. If you like hunting typos, open season begins when I announce this project on NA-DIGEST.
   2. Julia programming style, with the exception of correctness and performance. I know this is not fully idiomatic Julia, am working on it, and getting better. As I said in the introduction, I have traded a lot of abstraction for clarity. That means clairity for the novice. There may be something more abstract coming at the end of the project, but that is far away from now.
      1. I am also an old guy and the final product will reflect the Fortran __66__ I was raised on. That's show biz. 
           1.  Fortran + Julia = __Foolia__
   3. Organization of the repo. I'm still thinking this through. The important thing is that it make sense for the print book. I must do this work with the publisher.
   4. Questions like "Why isn't Trotsky's method in here?" If you object to an algorithmic choice, you'll have to be content to know that I have thought about the algorithm mix pretty carefully, have a clear vision for this project, and understand this field fairly well. 
   5. Questions like "Why doesn't SIAMFANLequations.jl look/work/smell like and/or use DasKapital.jl?" The reasons are that
      1. I am neither Karl nor Groucho,
      2. this project has a different mission, and 
      3. __I am working hard to limit depencencies__. 
   6. Philosophy, politics, opinions, invitations to debates, ...
 


## Core References and Documentation

The best documentation for this package will be the [notebook](https://github.com/ctkelley/NotebookSIAMFANL) and the print book. They will have detailed algorithmic descriptions, examples for you to play with, and guidance on tweaking the algorithmic paramenters to solve your problems. The notebook will be built in parallel with the print book and the content will be __roughly__ the same. The differences will be to accommodate the two formats. For example, docstrings need some work after the map from notebook to print and notebook has to make sense as an interactive resource.

I've also used [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) with this package. Click the badge
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

The solvers are designed to be stand-alone codes. The reason for this is the education mission of the project. I want the codes to be as easy to understand as possible. I have deliberately sacrificed a lot of abstraction and some performance in this effort. The reward for the reader (ie you) is that the algorithmic parameters are completely exposed so  you can play with them. At the end I hope to write a wrapper for all this that hides the parameters, but the stand-alone, keyword-infested codes are what you need if you want to really understand how these methods work. My students became experts in this field by fiddling with the Matlab version of these solvers.

The linear solvers are tuned to communicate well with nonlinear solvers. My old Matlab codes are a good illustration of this idea. My [new Mablab codes](https://ctk.math.ncsu.edu/knl.html) were designed in response to the need to do this better than I had been. In particular, the linear solver and the matrix-vector/preconditioner-vector product function need information on the nonlinear iteration and any precomputed data. While I could use global variables (and did in Kel95) and put these things in a module to simplify the interface, I won't do that anymore. Global varaibles make debugging harder and break parallelism. I like to avoid them. 

The algorithms, listed by book chapter will be

   - Chapter 1: Newton-Armijo and Pseudo-transient continuation for scalar equations: nsolsc.jl and ptcsolsc.jl
        - Codes: __Done!__, Notebook: __done!__
   - Chapter 2: Newton-Armijo and Pseudo-transient continuation for systems with direct linear solvers: nsold.jl and ptcd.jl
        - Codes: __Done!__, Notebook: __done!__
   - Chapter 3: Newton-Armijo and Pseudo-transient continuation for systems with iterative linear solvers: enable for nsoli.jl and ptcsoli.jl
       - nsoli.jl __done__
       - Linear solver(s): kl_gmres.jl (__90% done__) and kl_bicgstab.jl:
   - Chapter 4: Anderson acceleration: aasol.jl __Does Matlab code count as partially done?__
   - Chapter 5: Broyden's method: brsol.jl __0% done, but won't take long once I get started. I will do it the right way (ie from (Kel95)).__
   
   
## Test Problems

You'll need the TestProblems and examples submodules to run the notebook. To get those type 

```using SIAMFANLEquations.TestProblems```

and 

```using SIAMFANLEquations.Examples``` 

in the REPL or run the first code cell in the notebook 

```include("fanote_init.jl")```

There are two kinds of test problems. The ones you care about are the ones that I use in the print book and notebook to demonstrate the algorithms. The "inside baseball" problems are the ones I __only__ use for CI. They only appear in the /test directory. If you don't know or care what CI is, be happy.

## Citations 
Cite the package, print book and notebook like this. 
```
@misc{ctk:siamfanl,
title="{SIAMFANLEquations.jl}",
author="C. T. Kelley",
year=2021,
note="Julia Package",
doi="10.5281/zenodo.4284807",
url="https://github.com/ctkelley/SIAMFANLEquations.jl"
}

@misc{ctk:fajulia,
author="C. T. Kelley",
title="{Solving Nonlinear Equations with Iterative Methods:
Solvers and Examples in Julia}",
year=2021,
note="Unpublished book ms, under contract with SIAM"
}

@misc{ctk:notebooknl,
title="{Notebook for Solving Nonlinear Equations with Iterative Methods:
Solvers and Examples in Julia}",
author="C. T. Kelley",
year=2021,
note="IJulia Notebook",
url="https://github.com/ctkelley/NotebookSIAMFANL",
doi="10.5281/zenodo.4284687"
}
```

## Changes

### Updates since 0.2.3

- **0.3.0 is the current release.** It has Newton-GMRES (nsoli.jl) and GMRES (kl_gmres). 

- 0.3.1 happens when the notebook is closer to finished. As of 3/16/21
  - ptcsoli is working and covered by CI. A few examples are done.
  - restarted GMRES is working and in CI; need a few more examples for the notebook/real book
  - There's a lot of writing to be done before I can push 0.3.1
  
- Small things
   - Default side for preconditer is now __"right"__. See section 3.1.3 for the story on this.
   - Default forcing term is still constant __eta = .1__. This could change at any time and I've been careful to specify it completely in the examples.

### What's after 0.3.1?
 
 - 0.3.2 goes out when the writing is mostly done and BiCGSTAB is finished. 

 - 0.3.3 goes out when Chapter 3 is finished. I'm hoping for sometime in May. 

 - 0.4.0 is Anderson acceleration.

- 0.5.0 is Broyden's method

- 0.6.0 is the version I'll announce on NA-Digest, expect 0.5.x for x=0, ..., 9 before this happens.

   
## Funding

This project was partially supported by
1. Army Research Office grant W911NF-16-1-0504
2. National Science Foundation Grants
   1. OAC-1740309
   2. DMS-1745654
   3. DMS-1906446
3. Department of Energy grant DE-NA003967
   
Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author and
do not necessarily reflect the views of the National
Science Foundation, the Department of Energy,
or the Army Research Office.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ctkelley.github.io/SIAMFANLEquations.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ctkelley.github.io/SIAMFANLEquations.jl/stable

[build-status-img]: https://github.com/ctkelley/SIAMFANLEquations.jl/workflows/CI/badge.svg
[build-status-url]: https://github.com/ctkelley/SIAMFANLEquations.jl/actions

[codecov-img]: https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl



