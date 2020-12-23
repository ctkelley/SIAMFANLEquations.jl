module SIAMFANLEquations
using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SuiteSparse
using BandedMatrices
using Printf

export nsolsc
export ptcsolsc
export ptcsol
export nsol
export secant

include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
include("Tools/ManageStats.jl")
include("Tools/PTCTools.jl")
include("Chapter1/nsolsc.jl")
include("Chapter1/ptcsolsc.jl")
include("Chapter1/secant.jl")
include("Solvers/ptcsol.jl")
include("Solvers/nsol.jl")
include("PlotsTables/printhist.jl")
include("Tools/test_evaljac.jl")
include("Tools/NewtonIterationManagement.jl")

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
    beaminit,
    ptctest,
    pdeF!,
    pdeJ!,
    pdeinit,
    pdegminit,
    fishinit,
    fish2d,
    Pfish2d,
    Lap2d,
    Dx2d,
    Dy2d,
    solexact,
    l2dexact,
    dxexact,
    dyexact

include("TestProblems/Scalars/fcos.jl")
include("TestProblems/Scalars/fpatan.jl")
include("TestProblems/Scalars/spitchfork.jl")
include("TestProblems/Scalars/linatan.jl")
include("TestProblems/Scalars/ftanx.jl")
include("TestProblems/Systems/simple!.jl")
include("TestProblems/Systems/Fbvp!.jl")
include("TestProblems/Systems/FBeam!.jl")
include("TestProblems/Systems/Hequation.jl")
include("TestProblems/Systems/EllipticPDE.jl")
include("TestProblems/Systems/PDE_Tools.jl")
end

module Examples
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using LinearAlgebra
using BandedMatrices

export ptcBeam
export ivpBeam
export BVP_solve
export nsolheq
export NsolPDE

include("Examples/ptcBeam.jl")
include("Examples/ivpBeam.jl")
include("Examples/BVP_solve.jl")
include("Examples/NsolPDE.jl")
include("Examples/Internal/nsolheq.jl")
end

end # module
