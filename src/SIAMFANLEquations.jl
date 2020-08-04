module SIAMFANLEquations
using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SuiteSparse
using Printf

export nsolsc
export parab3p
export ptcsolsc
export printhist
export nsold
#export difffp
#export fpeval_newton
#export armijosc

include("Tools/parab3p.jl")
include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
include("Tools/FunctionDerivativeEvals.jl")
include("Tools/ManageStats.jl")
include("Chapter1/nsolsc.jl")
include("Chapter1/ptcsolsc.jl")
include("Chapter2/nsold.jl")
include("PlotsTables/printhist.jl")
include("Tools/test_evaljac.jl")

module TestProblems
using LinearAlgebra
using SparseArrays
using SuiteSparse
using AbstractFFTs
using FFTW
using Printf

export
#Functions
fcos,
fpatan,
spitchfork,
linatan,
sptestp,
sptest,
ftanx,
ftanxp,
heqinit,
heqf!,
heqJ!,
simple!,
jsimple!,
heqbos!,
setc!,
chandprint


include("TestProblems/Scalars/fcos.jl")
include("TestProblems/Scalars/fpatan.jl")
include("TestProblems/Scalars/spitchfork.jl")
include("TestProblems/Scalars/linatan.jl")
include("TestProblems/Scalars/ftanx.jl")
include("TestProblems/Chapter2/simple!.jl")
include("TestProblems/Heq4nsold.jl")

using .Heq4nsold
end


end # module
