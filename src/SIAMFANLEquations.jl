module SIAMFANLEquations
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Printf

export nsolsc
export parab3p
export difffp
export fpeval_newton
export ptcsolsc
export printhist
export armijosc

#mutable struct ItStats
#ifun::Array{Int64,1}
#ijac::Array{Int64,1}
#iarm::Array{Int64,1}
#history::Array{T,1} where T<:Real
#end

include("Tools/difffp.jl")
include("Tools/fpeval_newton.jl")
include("Tools/parab3p.jl")
include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionDerivativeEvals.jl")
include("Tools/ManageStats.jl")
include("Chapter1/nsolsc.jl")
include("Chapter1/ptcsolsc.jl")
include("PlotsTables/printhist.jl")

module TestProblems
using LinearAlgebra
using SparseArrays
using SuiteSparse
using AbstractFFTs
using FFTW

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

include("TestProblems/Scalars/fcos.jl")
include("TestProblems/Scalars/fpatan.jl")
include("TestProblems/Scalars/spitchfork.jl")
include("TestProblems/Scalars/linatan.jl")
include("TestProblems/Scalars/ftanx.jl")
end


end # module
