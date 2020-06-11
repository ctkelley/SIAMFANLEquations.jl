module SIAMFANLEquations
using LinearAlgebra
using SparseArrays
using SuiteSparse
using Printf

export nsolsc
export parab3p
export difffp
export fpeval_newton
export ptcsc
export printhist

include("Tools/difffp.jl")
include("Tools/fpeval_newton.jl")
include("Tools/parab3p.jl")
include("Chapter1/nsolsc.jl")
include("Chapter1/ptcsc.jl")
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
