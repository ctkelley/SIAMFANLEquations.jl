| **Documentation**                                                               | **Build Status**                                                                                | **DOI**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------- |
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][build-status-img]][build-status-url] [![][codecov-img]][codecov-url] | [![DOI](https://zenodo.org/badge/256312455.svg)](https://zenodo.org/badge/latestdoi/256312455) |

[![SIAMFANLEquaitons Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/SIAMFANLEquations)](https://pkgs.genieframework.com?packages=SIAMFANLEquations)

# SIAMFANLEquations 


## The archival version 1.0 from the date of publication is in the [FA20 branch](https://github.com/ctkelley/SIAMFANLEquations.jl/tree/FA20).


## This is the Julia package for my shiny new orange book

<img width = 400, src="https://user-images.githubusercontent.com/10243067/184647769-d9d51ee9-79f0-48ba-96a4-b9ed2a66cdfa.jpg">

# [Solving Nonlinear Equations with Iterative Methods: <br> Solvers and Examples in Julia](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635)

## [C. T. Kelley](https://ctk.math.ncsu.edu)

The book is finished and this project is __DONE__. So I take the sacred book author oath ...
  - I will only make updates to the package and notebooks to fix bugs or typos. 
  - I will not be adding new functionality to this package or new material to the notebooks. 
  - I will make no changes to the user interface for the codes in the package.

This is a sequel to my book 

(Kel03) C. T. Kelley, [***Solving Nonlinear Equations with Iterative Methods:***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM, Philadelphia, 2003.

Hence the notebook and this package all have SIAMFANL in their names.

The new book has a different algorithm mix and the solvers and examples are in Juila. The project has three parts.

   1. A print book: (Kel22) C. T. Kelley, [***Solving Nonlinear Equations with Newton's Method: Solvers and Examples in Julia***](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635), Fundamentals of Algorithms 20, SIAM, Philadelphia, 2022.    
   3. [A suite of IJulia notebooks](https://github.com/ctkelley/NotebookSIAMFANL) (open source, MIT License, Creative Commons License)
      The latest releases of the notebook suite and package run correctly. The notebooks and package from the master branches also run correctly
      together. Bug fixes prior to 1.0 may, with an absurdly low probablilty, break things in older releases. 
   3. This package (MIT License)<br>

Content changes from (Kel03):

- New solvers: __pseudo-transient continuation__ and __Anderson acceleration__
- Deletions: __Broyden's method__ 
    - Quasi-Newton methods are not used much for nonlinear equations any more. Newton-Krylov has taken over.
- New Case Studies chapter
   
## Readme Contents:
- [Mission](#package-Mission)
- [Installation](#installation)
- [Reporting bugs: __Please__ No Pull Requests](#pull-Requests)
- [Core References and Documentation](#Core-References-and-Documentation)
- [Algorithms and Solvers](#Algorithms-and-Solvers)
- [About the test problems](#Test-Problems-and-the-notebook)
- [How to cite this stuff](#Citations)
- [Book FAQs](#FAQs)
- [Funding](#Funding) 

## Package Mission

This package is designed and built to support a book project. So the solvers and examples reinforce the algorithmic discussion in the book. General purpose packages have a different mission. 

## Installation: 


- Your best bet is to __use the latest version of Julia__  (currently 1.9.2) with the notebooks and the package.
- If you must use old stuff, use LTS 1.6.7 and up with this thing!!!
- Please do not use any non-LTS version earlier than 1.8. The notebook kernel is now 1.9.

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
to run the test problems. Then there are the examples you get with
```
using SIAMFANLEquations.Examples
```
for the unit tests, the examples in the book, and the notebook.


## Pull Requests

My favorite thing about book projects is that they are not open-ended. They get finished. For example, take 
[this book](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635) ... please.

__Even after publication, I like bug reports; I need bug reports__, but ...

__Please, please__, do not send me PRs. If you find 
   1.  a bug (programming or performance) in the codes,
   2.  errors and/or typos in the notebooks/docstrings/readme
   3. confusion, lack of clarity, or __errors in the installation instructions__,
       1. I would __really like__ some Windows users to try this stuff, especially the notebooks.
   4. something I could do in the codes to help you do your work ...
       1. that won't break other stuff, which includes the connection between the book and the package,
       2. or eat up lots of time,
  
 Please  ... 

- tell me the old fashioned way with email to tim_kelley@ncsu.edu 
- or open an issue.

This is a book project and I need to put all changes in by hand so I'll have muscle memory about what's going on. If there is a second printing I can fix things in the print/pdf books and will fix things in real time (more or less) in the codes and notebooks.

I have limited bandwidth, __so please do not send me email or open issues about__ ...

   1. Julia programming style, with the exception of correctness and performance. I know this is not fully idiomatic Julia. I got somewhat better as the project progressed. As I said in the introduction, I have traded a lot of abstraction for clarity. That means clairity for the novice. 
      1. I am also an old guy and the final product will reflect the Fortran __66__ I was raised on. That's show biz. 
           1.  Fortran + Julia = __Foolia__
   3. Questions like "Why isn't Trotsky's method in here?" If you object to an algorithmic choice, you'll have to be content to know that I thought about the algorithm mix pretty carefully, had a clear vision for this project, and understand this field fairly well. 
   4. Questions like "Why doesn't SIAMFANLEquations.jl look/work/smell like and/or use DasKapital.jl?" The reasons are that
      1. I am neither Karl nor Groucho,
      2. this project has a different mission, and 
      3. __I worked hard to limit depencencies__. 
   5. Philosophy, politics, opinions, invitations to debates, ...
   6. Organization of the repo, names of functions, API, or anything else that is now __frozen for the book__.


## Core References and Documentation

The best documentation for this package lives in the [notebook](https://github.com/ctkelley/NotebookSIAMFANL) and the print book. They have detailed algorithmic descriptions, examples for you to play with, and guidance on tweaking the algorithmic paramenters to solve your problems. The notebook was built in parallel with the print book and the content is __roughly__ the same. The differences are mostly to accommodate the two formats. For example, docstrings need some work after the map from notebook to print and the notebook has to make sense as an interactive resource.

I've also used [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) with this package. Click the badge
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/stable)
to get the documentation from the latest release. The documenter files have the headers for the solvers and some of the test problems. I continue to work on the docs and they will get better, but will never be as good as the notebook.

This book will not cover theory in detail (ie no proofs). My two books on nonlinear equations

(Kel95) C. T. Kelley, [***Iterative Methods for Linear and Nonlinear Equations***](https://my.siam.org/Store/Product/viewproduct/?ProductId=862) , Frontiers in Applied Mathematics 16,  SIAM, Philadelphia, 1995

and

(Kel03) C. T. Kelley, [***Solving Nonlinear Equations with Newton's Method***](https://my.siam.org/Store/Product/viewproduct/?ProductId=841) , Fundamentals of Algorithms 1, SIAM, Philadelphia, 2003

describe the classic Newton and Newton-Krylov algorithms. Kel95 has the theory. This project is a sequel to Kel03. Kel03 is Matlab-centric
and will remain in print.

A recent Acta Numerica paper has everything

(Kel18) C. T. Kelley, ***Numerical Methods for Nonlinear Equations***, Acta Numerica 27 (2018), pp 207--287. https://doi.org/10.1017/S0962492917000113

The references I use for theory of pseudo-transient continuation and Anderson acceleration are

(KK98) C. T. Kelley and D. E. Keyes, ***Convergence Analysis of Pseudo-Transient Continuation***, SIAM Journal on Numerical Analysis 35 (1998), pp 508-523. https://doi.org/10.1137/S0036142996304796

(TK15) A. Toth and C. T. Kelley, ***Convergence Analysis for Anderson Acceleration***, SIAM Journal on Numerical Analysis 53, (2015), pp 805-819. https://doi.org/10.1137/130919398

## Algorithms and Solvers

The solvers are designed to be stand-alone codes. The reason for this is the education mission of the project. I want the codes to be as easy to understand as possible. I have deliberately sacrificed a lot of abstraction and some performance in this effort. The reward for the reader (ie you) is that the algorithmic parameters are completely exposed so  you can play with them. Someday, not soon, I may write a wrapper for all this that hides the parameters as a separate package. However, the stand-alone, keyword-infested codes are what you need if you want to really understand how these methods work. My students became experts in this field by fiddling with the Matlab version of these solvers.

The linear solvers are tuned to communicate well with nonlinear solvers. My old Matlab codes are a good illustration of this idea. My [new Mablab codes](https://ctk.math.ncsu.edu/knl.html) were designed in response to the need to do this better than I had been. In particular, the linear solver and the matrix-vector/preconditioner-vector product function need information on the nonlinear iteration and any precomputed data. While I could use global variables (and did in Kel95) and put these things in a module to simplify the interface, I won't do that anymore. Global varaibles make debugging harder and break parallelism. I like to avoid them. 

The algorithms, listed by book chapter are

   - Chapter 1: Newton-Armijo and Pseudo-transient continuation for scalar equations: __nsolsc.jl__ and __ptcsolsc.jl__
       
   - Chapter 2: Newton-Armijo and Pseudo-transient continuation for systems with direct linear solvers: __nsol.jl__ and __ptcsol.jl__
       
   - Chapter 3: Newton-Armijo and Pseudo-transient continuation for systems with iterative linear solvers: __nsoli.jl__ and __ptcsoli.jl__
       
   - Chapter 4: Anderson acceleration: __aasol.jl__ 
        
   - Chapter 5: Case studies:  __Conductive-Radiative heat transfer__ and __Continuation for H-equation.__ 
   
   
## Test Problems and the notebook

You'll need the TestProblems and Examples submodules to run the notebook. To get those type 

```using SIAMFANLEquations.TestProblems```

and 

```using SIAMFANLEquations.Examples``` 

in the REPL or run the first code cell in the notebook 

```include("fanote_init.jl")```

There are two kinds of test problems. The ones you care about are the ones that I use in the print book and notebook to demonstrate the algorithms. The "inside baseball" problems are the ones I __only__ use for CI. They only appear in the /test directory. If you don't know or care about what CI is, be happy.

## Citations 
Cite the package, print book and notebook like this. 
```
@misc{ctk:siamfanl,
title="{SIAMFANLEquations.jl}",
author="C. T. Kelley",
year=2022,
note="Julia Package",
doi="10.5281/zenodo.4284807",
url="https://github.com/ctkelley/SIAMFANLEquations.jl"
}

@book{ctk:fajulia,
author="C. T. Kelley",
title="{Solving Nonlinear Equations with Iterative Methods:
Solvers and Examples in Julia}",
year=2022,
publisher="SIAM",
address="Philadelphia",
series="Fundamentals of Algorithms",
number=20
}

@misc{ctk:notebooknl,
title="{Notebook for Solving Nonlinear Equations with Iterative Methods:
Solvers and Examples in Julia}",
author="C. T. Kelley",
year=2022,
note="IJulia Notebook",
url="https://github.com/ctkelley/NotebookSIAMFANL",
doi="10.5281/zenodo.4284687"
}
```

## FAQs
1. What kind of book is this?

    - It's an orange book.

2. What is this book about?

     - It's about 200 pages.

3. Have you written any other amazing books?

     - [Yes.](https://ctk.math.ncsu.edu/lv/books.html)

   
## Funding

This project was partially supported by

1. National Science Foundation Grants
   1. OAC-1740309
   2. DMS-1745654
   3. DMS-1906446
2. Department of Energy grant DE-NA003967
3. Army Research Office grant W911NF-16-1-0504
   
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



