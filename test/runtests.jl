using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test
using LinearAlgebra
include("Chapter1/nsolsc_solution_test.jl")
include("Chapter1/ptcsc_test.jl")
@test nsolsc_solution_test()
@test ptcsc_test()
