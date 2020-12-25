using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using SIAMFANLEquations.Examples
using Test
using LinearAlgebra

include("Chapter1/nsolsc_solution_test.jl")
include("Chapter1/ptcsolsc_test.jl")
include("Chapter2/basic2d_test.jl")
include("Chapter2/heq_lu_test.jl")
include("Chapter2/bvp_test.jl")
include("Chapter2/beam_test.jl")
include("Chapter2/pde_lin_test.jl")
include("Chapter2/nsolpde_test.jl")
include("Chapter3/gmres_test.jl")
@testset "Scalar Equations: Chapter 1" begin
    @test nsolsc_solution_test()
    @test ptcsolsc_test()
end
@testset "nsol and ptcsol: Chapter 2" begin
    @test basic2d_test()
    @test bvp_test(201)
    @test beam_test()
    @test heq_lu_test()
    @test pde_lin_test(31)
    @test nsolpde_test(31)
end
@testset "Krylov solvers: Chapter 3" begin
    @test gmres_test()
end
