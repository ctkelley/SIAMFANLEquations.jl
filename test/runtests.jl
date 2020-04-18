using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test
using LinearAlgebra
include("nsolsc_solution_test.jl")
include("ptcsc_test.jl")
@test nsolsc_solution_test()
@test ptcsc_test()
