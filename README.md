[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ctkelley.github.io/SIAMFANLEquations.jl/dev)
[![Build Status](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl.svg?branch=master)](https://travis-ci.com/ctkelley/SIAMFANLEquations.jl)
[![codecov](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ctkelley/SIAMFANLEquations.jl)
# SIAMFANLEquations
This the start to the package for the book. It is not ready yet so ...

Please stop reading and go away.

If you're me, here are the notes

1. The github address for the package is https://github.com/ctkelley/SIAMFANLEquations.jl
   1. The way to install it is to type
   
      add https://github.com/ctkelley/SIAMFANLEquations.jl
      
      at the pkg prompt.
      
   2. The way to make changes and get them back on github is 
         1. Start Julia in the package directory with 
         
            Julia --project==.
            
            This makes "using SIAMFANLEquations" point to the directory rather than the version in .julia/packages
         2. Make the changes and push them to the repo
         3. quit Julia
         4. restart Julia and update the package with pkg
         
