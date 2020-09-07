module SIAMFANLEquations
using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SuiteSparse
using BandedMatrices
using Printf

export nsolsc
export ptcsolsc
export nsold

include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
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
using BandedMatrices
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
    chandprint,
    bvpinit,
    Fbvp!,
    Jbvp!,
    FBeam!,
    FBeamtd!,
    BeamJ!,
    BeamtdJ!,
    beaminit
    


include("TestProblems/Scalars/fcos.jl")
include("TestProblems/Scalars/fpatan.jl")
include("TestProblems/Scalars/spitchfork.jl")
include("TestProblems/Scalars/linatan.jl")
include("TestProblems/Scalars/ftanx.jl")
include("TestProblems/Systems/simple!.jl")
include("TestProblems/Systems/Fbvp!.jl")
include("TestProblems/Systems/FBeam!.jl")
include("TestProblems/Systems/Heq4nsold.jl")

end


end # module
