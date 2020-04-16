using Documenter, ScalarEquations, DocumenterLaTeX, DocumenterTools
push!(LOAD_PATH,"../src/")
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end
makedocs(
   sitename="ScalarEquations.jl",
   authors="C. T. Kelley",
   format = Documenter.HTML(
   prettyurls = !("local" in ARGS),
   ),
   pages = Any[
     "Home" => "index.md",
     "Scalar Equations" => Any[
        "Scalar.md"
      ],
     "Scalar Equations Functions" => Any[
       "functions/functions.md",
       "functions/atan_test.md",
      ]
    ]
)
#mkpath("build/images")
#cp("../images/newton_big.jpg","build/images/newton_big.jpg",force=true)
deploydocs(
     repo="github.com/ctkelley/ScalarEquations.jl.git",
     target="build",
     devurl="stable"
)

