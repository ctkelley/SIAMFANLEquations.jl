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
export ptcsoli
export nsol
export nsoli
export nofact
export secant
export armijosc
export kl_gmres
export kl_bicgstab
export Katv
export kstore
export Orthogonalize!
export EvalF!
export solhistinit

include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
include("Tools/ManageStats.jl")
include("Tools/IterationInit.jl")
include("Tools/ErrorTest.jl")
include("Tools/NewtonKrylov_Tools.jl")
include("Tools/PTCTools.jl")
include("Tools/PTCToolsi.jl")
include("Solvers/Chapter1/nsolsc.jl")
include("Solvers/Chapter1/ptcsolsc.jl")
include("Solvers/Chapter1/secant.jl")
include("Solvers/ptcsol.jl")
include("Solvers/ptcsoli.jl")
include("Solvers/nsol.jl")
include("Solvers/nsoli.jl")
include("Solvers/LinearSolvers/kl_gmres.jl")
include("Solvers/LinearSolvers/kl_bicgstab.jl")
include("Solvers/LinearSolvers/Orthogonalize!.jl")
include("PlotsTables/printhist.jl")

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
    JVsimple,
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
    Jvec2d,
    pdeinit,
    pdegminit,
    fishinit,
    fish2d,
    sintv,
    isintv,
    Pfish2d,
    Pvec2d,
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

export ptciBeam
export ptcBeam
export ivpBeam
export BVP_solve
export nsolheq
export NsolPDE
export NsoliPDE

include("Examples/ptciBeam.jl")
include("Examples/ptcBeam.jl")
include("Examples/ivpBeam.jl")
include("Examples/BVP_solve.jl")
include("Examples/NsolPDE.jl")
include("Examples/NsoliPDE.jl")
include("Examples/Internal/nsolheq.jl")
end

end # module
