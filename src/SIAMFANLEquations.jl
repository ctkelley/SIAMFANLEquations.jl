module SIAMFANLEquations

#using PyPlot
#using LaTeXStrings

export nsolsc
export parab3p
export difffp
export fpeval_newton
export ptcsc

include("Tools/difffp.jl")
include("Tools/fpeval_newton.jl")
include("Tools/parab3p.jl")
include("nsolsc.jl")
include("ptcsc.jl")

module TestProblems
using LinearAlgebra
using SparseArrays
using SuiteSparse

export
#Functions
fcos,
fpatan,
spitchfork,
linatan,
sptestp,
sptest,
ftanx,
ftanxp

include("TestProblems/fcos.jl")
include("TestProblems/fpatan.jl")
include("TestProblems/spitchfork.jl")
include("TestProblems/linatan.jl")
include("TestProblems/ftanx.jl")
end


end # module
