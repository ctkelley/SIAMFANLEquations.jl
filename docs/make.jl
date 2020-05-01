using Documenter, SIAMFANLEquations, DocumenterLaTeX, DocumenterTools
push!(LOAD_PATH,"../src/")
makedocs(sitename="SIAMFANLEquations.jl",
authors="C. T. Kelley",
format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
pages = Any[
     "Home" => "index.md",
     "Scalar Equation Solvers" => Any[
       "functions/nsolsc.md",
       "functions/ptcsc.md",
      ]
]
)
deploydocs(
     repo="github.com/ctkelley/SIAMFANLEquations.jl/tree/dev"
)
