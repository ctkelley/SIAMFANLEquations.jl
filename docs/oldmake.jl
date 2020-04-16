using Documenter, ScalarEquations, DocumenterLaTeX, DocumenterTools
push!(LOAD_PATH,"../src/")
makedocs(sitename="ScalarEquations.jl",
authors="C. T. Kelley",
format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
pages = Any[
     "Home" => "index.md",
     "Scalar Equations" => Any[
        "Scalar.md"
      ],
     "Scalar Equations Functions" => Any[
       "functions/functions.md",
       "functions/sptc.md",
      ]
]
)
deploydocs(
     repo="github.com/ctkelley/ScalarEquations.jl.git"
)
