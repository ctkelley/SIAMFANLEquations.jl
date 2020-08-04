using Documenter, SIAMFANLEquations, DocumenterLaTeX, DocumenterTools
push!(LOAD_PATH,"../src/")
makedocs(sitename="SIAMFANLEquations.jl",
authors="C. T. Kelley",
format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
pages = Any[
     "Home" => "index.md",
     "Scalar Equations" => Any[
       "functions/nsolsc.md",
       "functions/ptcsolsc.md",
       ],
     "Systems with Direct Linear Solvers" => Any[
       "functions/nsold.md",
      ]
     
]
)
deploydocs(
     repo="github.com/ctkelley/SIAMFANLEquations.jl.git"
)
