using Test
using LinearAlgebra

using HelmholtzLAP

@testset "HelmholtzLAP" begin
    @testset "module interface" begin
        @test isdefined(HelmholtzLAP, :assemble_fem)
        @test :assemble_fem in names(HelmholtzLAP)
    end

    @testset "fem" begin
        meshfile = joinpath(@__DIR__, "..", "circle", "geo-circle.msh")
        @test isfile(meshfile)

        fem = HelmholtzLAP.assemble_fem(meshfile)

        @test propertynames(fem) == (:model, :V0, :U, :Ω, :dΩ, :K, :M)

        K = fem.K
        M = fem.M

        # The stiffness and mass matrices must be nonempty square matrices
        # with matching dimensions.
        @test size(K, 1) > 0
        @test size(K, 1) == size(K, 2)
        @test size(M) == size(K)

        # Symmetric bilinear forms should produce numerically symmetric matrices.
        @test isapprox(K, K'; rtol=1e-12, atol=1e-12)
        @test isapprox(M, M'; rtol=1e-12, atol=1e-12)

        # Use one nonzero vector for a basic check of positive quadratic forms.
        x = ones(size(K, 1))
        @test dot(x, K * x) > 0
        @test dot(x, M * x) > 0

        eigenpairs = HelmholtzLAP.compute_eigenpairs(K, M; nev=4, which=:SM)

        @test propertynames(eigenpairs) == (:λ, :Φ)

        λ = eigenpairs.λ
        Φ = eigenpairs.Φ

        @test length(λ) == 4
        @test size(Φ, 1) == size(K, 1)

        HelmholtzLAP.show_eigens(λ)
    end

end
