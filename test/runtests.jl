using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test
using LinearAlgebra
include("Chapter1/nsolsc_solution_test.jl")
include("Chapter1/ptcsolsc_test.jl")
@testset "Scalar Equations: Chapter 1" begin
   @test nsolsc_solution_test()
   @test ptcsolsc_test()
end
