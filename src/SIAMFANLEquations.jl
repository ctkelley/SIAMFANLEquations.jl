module SIAMFANLEquations

#using PyPlot
#using LaTeXStrings

export nsolsc
export difffp
export fpeval_newton
export fpatan
export fcos
export ptcsc
export sptestp
export sptest
export linatan

include("Tools/difffp.jl")
include("Tools/fpeval_newton.jl")
include("nsolsc.jl")
include("ptcsc.jl")
include("TestProblems/fcos.jl")
include("TestProblems/fpatan.jl")
include("TestProblems/spitchfork.jl")
include("TestProblems/linatan.jl")


end # module
