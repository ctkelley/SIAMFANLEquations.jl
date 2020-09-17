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
export ptcsold
export ptcsol
export PTCinit

include("Tools/armijo.jl")
include("Tools/PrintError.jl")
include("Tools/FunctionJacobianEvals.jl")
include("Tools/ManageStats.jl")
include("Solvers/PTCTools.jl")
include("Chapter1/nsolsc.jl")
include("Chapter1/ptcsolsc.jl")
include("Chapter2/nsold.jl")
include("Solvers/ptcsol.jl")
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
include("TestProblems/Systems/Hequation.jl")


end


module NotebookSIAMFANL

using PyPlot
using LaTeXStrings
using Printf
using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems

export atan_test
export ftan_test
export ptc_scalar_example
export tab1dot1
export fig1dot1
export fig1dot2
export fig1dot3
export fig1dot3b
export fig1dot4
export fig1dot5
export fig1dot6
export fig1dot7
export fig2dot1
export MyFunction
export PitchFork1

export plothist

include("Notebook/Chapter1/fig1dot1.jl")
include("Notebook/Chapter1/fig1dot2.jl")
include("Notebook/Chapter1/fig1dot3.jl")
include("Notebook/Chapter1/fig1dot3b.jl")
include("Notebook/Chapter1/fig1dot4.jl")
include("Notebook/Chapter1/fig1dot5.jl")
include("Notebook/Chapter1/fig1dot6.jl")
include("Notebook/Chapter1/fig1dot7.jl")
include("Notebook/Chapter1/tab1dot1.jl")
include("Notebook/Chapter1/MyFunction.jl")
include("Notebook/Chapter1/PitchFork1.jl")
include("Notebook/Chapter2/fig2dot1.jl")

include("Notebook/Tools/plothist.jl")

end # module NotebookSIAMFANL

end # module
