module SIAMFANLEquations
using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SuiteSparse
using BandedMatrices
using LaTeXStrings
#using QuadGK
using Printf

# Export the nonlinear solvers
export nsolsc
export ptcsolsc
export ptcsol
export ptcsoli
export nsol
export nsoli
export aasol
export secant
# nofact is the factorization that does nothing. It is
# a legal option for nsol and ptcsol and I must export it.
export nofact
# Export the linear solvers
export kl_gmres
export kl_bicgstab
# A couple functions the solvers need to manage storage. These are
# from src/Tools/NewtonKrylov_Tools.jl
#
# kstore gets the vectors GMRES needs internally and makes room to
# copy the initial iterate and right side. I use this in the heat
# transfer problem in Chapter 5.
#
export kstore
#
# knl_init preallocates the vectores nsoli and ptcsoli use internally.
# I need to export this for the continuation code in Chapter 5.
#
export nkl_init
#
#
#export knlstore
#export EvalF!
#export solhistinit
#export armijosc
#export Katv
#export Orthogonalize!

include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
include("Tools/ManageStats.jl")
include("Tools/IterationInit.jl")
include("Tools/ErrorTest.jl")
include("Tools/NewtonKrylov_Tools.jl")
include("Tools/PTCTools.jl")
include("Tools/AA_Tools.jl")
include("Solvers/Chapter1/nsolsc.jl")
include("Solvers/Chapter1/ptcsolsc.jl")
include("Solvers/Chapter1/secant.jl")
include("Solvers/ptcsol.jl")
include("Solvers/ptcsoli.jl")
include("Solvers/nsol.jl")
include("Solvers/nsoli.jl")
include("Solvers/aasol.jl")
include("Solvers/LinearSolvers/kl_gmres.jl")
include("Solvers/LinearSolvers/kl_bicgstab.jl")
include("Solvers/LinearSolvers/Orthogonalize!.jl")
#include("PlotsTables/printhist.jl")

module TestProblems
using SIAMFANLEquations
using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SuiteSparse
using BandedMatrices
using AbstractFFTs
using FFTW
using Printf
using QuadGK

export
    #Functions
    #    fcos,
    #    fpatan,
    spitchfork,
    #    linatan,
    sptestp,
    sptest,
    #    ftanx,
    #    ftanxp,
    heqinit,
    heqf!,
    heqJ!,
    HeqFix!,
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
    Lap1d,
    Dx2d,
    Dy2d,
    solexact,
    l2dexact,
    dxexact,
    dyexact,
    hardleft!,
    hardleftFix!,
    heat_init,
    sn_init,
    heat_fixed!,
    FCR_heat!,
    getrhs,
    AxB,
    transport_sweep!,
    heq_continue,
    knl_continue



include("TestProblems/Scalars/spitchfork.jl")
include("TestProblems/Systems/simple!.jl")
include("TestProblems/Systems/Fbvp!.jl")
include("TestProblems/Systems/FBeam!.jl")
include("TestProblems/Systems/Hequation.jl")
include("TestProblems/Systems/EllipticPDE.jl")
include("TestProblems/Systems/PDE_Tools.jl")
include("TestProblems/CaseStudies/CR_Heat.jl")
include("TestProblems/CaseStudies/knl_continue.jl")
include("TestProblems/CaseStudies/heq_continue.jl")
end

module Examples
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using LinearAlgebra
using LinearAlgebra.BLAS
using BandedMatrices

export ptciBeam
export ptcBeam
export ivpBeam
export BVP_solve
export nsolheq
export NsolPDE
export NsoliPDE
export PDE_aa

include("Examples/ptciBeam.jl")
include("Examples/ptcBeam.jl")
include("Examples/ivpBeam.jl")
include("Examples/BVP_solve.jl")
include("Examples/NsolPDE.jl")
include("Examples/NsoliPDE.jl")
include("Examples/PDE_aa.jl")
include("Examples/Internal/nsolheq.jl")
end

end # module
