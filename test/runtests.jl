using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using SIAMFANLEquations.Examples
using Test
using LinearAlgebra: LinearAlgebra, Diagonal, I, diagm, mul!, norm, qr, qr!
#using LinearAlgebra.BLAS
import SIAMFANLEquations.Orthogonalize!

include("Chapter1/nsolsc_solution_test.jl")
include("Chapter1/ptcsolsc_test.jl")
include("Chapter2/basic2d_test.jl")
include("Chapter2/heq_lu_test.jl")
include("Chapter2/bvp_test.jl")
include("Chapter2/beam_test.jl")
include("Chapter2/pde_lin_test.jl")
include("Chapter2/nsolpde_test.jl")
include("Chapter2/knowsdt_test.jl")
include("Chapter3/gmres_test.jl")
include("Chapter3/mgs_test.jl")
include("Chapter3/bicgstab_test.jl")
include("Chapter3/Krylov_pde_test.jl")
include("Chapter3/ptcKrylovTest.jl")
include("Chapter3/ptcKrylovTestB.jl")
include("Chapter3/nk_test.jl")
include("Chapter3/nk_pde.jl")
include("Chapter3/nk_heq.jl")
include("Chapter4/reldiff.jl")
include("Chapter4/alex_test.jl")
include("Chapter4/ci_pde_aa.jl")
include("Chapter4/heq_aa.jl")
include("Chapter4/linear_aa.jl")
include("Chapter5/transport_test.jl")
include("Chapter5/heat_test.jl")
include("Chapter5/heat_test2.jl")
include("Chapter5/continue_test.jl")


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
    @test knowsdt_test()
end
@testset "Newton-Krylov solvers: Chapter 3" begin
    @test nk_test()
    @test nk_pde()
    @test nk_heq()
    @test ptcKrylovTest()
    @test ptcKrylovTestB()
end
@testset "Krylov solvers: Chapter 3" begin
    @test gmres_test()
    @test mgs_test()
    @test bicgstab_test()
    @test gmres_test_pde(31)
    @test bicgstab_test_pde(31)
end
@testset "Anderson Acceleration: Chapter 4" begin
    @test ci_pde_aa()
    @test heq_aa()
    @test linear_aa()
    @test alex_test()
end
@testset "Case Studies: Chapter 5" begin
    @test transport_test()
    @test heat_test()
    @test heat_test2()
    @test continue_test()
end
