using Test
using LinearAlgebra

using HelmholtzLAP

@testset "HelmholtzLAP" begin
    @testset "module interface" begin
        public_functions = (
            :assemble_fem,
            :assemble_source_loads,
            :compute_eigenpairs,
            :show_eigenvalues,
            :find_resonant_indices,
            :check_source_orthogonality,
            :project_load_to_range,
            :solve_lap,
        )

        internal_functions = (
            :prepare_source_extension,
            :compute_particular_solution,
            :compute_lap_solution,
        )

        for function_name in public_functions
            @test isdefined(HelmholtzLAP, function_name)
            @test function_name in names(HelmholtzLAP)
        end

        for function_name in internal_functions
            @test isdefined(HelmholtzLAP, function_name)
            @test !(function_name in names(HelmholtzLAP))
        end
    end

    @testset "fem" begin
        meshfile = joinpath(@__DIR__, "test_circle", "geo-circle.msh")
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

        eigenpairs = HelmholtzLAP.compute_eigenpairs(fem; nev=4, which=:SM)

        @test propertynames(eigenpairs) == (:λ, :Φ)

        λ = eigenpairs.λ
        Φ = eigenpairs.Φ

        @test length(λ) == 4
        @test size(Φ, 1) == size(K, 1)

        HelmholtzLAP.show_eigenvalues(eigenpairs)

        resonant_indices = HelmholtzLAP.find_resonant_indices(eigenpairs, 2)
        @test resonant_indices == 2:3
        @test_throws ArgumentError HelmholtzLAP.find_resonant_indices(eigenpairs, 0)

        source(x, k) = k * x[1]
        source_k(x, k) = x[1]
        source_loads = HelmholtzLAP.assemble_source_loads(fem, source, source_k, 1.0)

        @test length(source_loads.b) == size(K, 1)
        @test length(source_loads.b_der) == size(K, 1)
    end

    @testset "LAP bundled inputs" begin
        fem = (K=[1.0 0.0; 0.0 2.0], M=[1.0 0.0; 0.0 1.0])
        Φres = reshape([1.0, 0.0], 2, 1)

        orthogonal_loads = (b=[0.0, 2.0], b_der=[4.0, 3.0])
        nonorthogonal_loads = (b=[5.0, 2.0], b_der=[4.0, 3.0])

        check = HelmholtzLAP.check_source_orthogonality(orthogonal_loads, Φres)
        @test check.is_orthogonal

        projection = HelmholtzLAP.project_load_to_range(fem, nonorthogonal_loads.b, Φres)
        @test projection.coefficients ≈ [5.0]
        @test projection.projected ≈ [0.0, 2.0]
        @test norm(Φres' * projection.projected) ≤ 1e-12

        frozen = HelmholtzLAP.prepare_source_extension(
            fem,
            nonorthogonal_loads,
            Φres,
            :frozen_correction,
        )
        @test frozen.b_effective ≈ [0.0, 2.0]
        @test frozen.b_der_effective == nonorthogonal_loads.b_der

        projected = HelmholtzLAP.prepare_source_extension(
            fem,
            nonorthogonal_loads,
            Φres,
            :projected_correction,
        )
        @test projected.b_effective ≈ [0.0, 2.0]
        @test projected.b_der_effective ≈ [0.0, 3.0]

        particular = HelmholtzLAP.compute_particular_solution(
            fem,
            orthogonal_loads,
            1.0,
            Φres;
            source_mode=:original,
            restart=2,
            maxiter=10,
            reltol=1e-12,
        )
        @test particular.u_special ≈ [0.0, 2.0]
        @test particular.linear_residual ≤ 1e-12

        lap = HelmholtzLAP.compute_lap_solution(fem, particular, Φres, 1.0)
        @test lap.coefficients ≈ [-2.0]
        @test lap.u_lap ≈ [-2.0, 2.0]
        @test lap.constraint_residual ≤ 1e-12

        solution = HelmholtzLAP.solve_lap(
            fem,
            orthogonal_loads,
            Φres,
            1.0;
            source_mode=:original,
            restart=2,
            maxiter=10,
            reltol=1e-12,
        )
        @test solution.particular.u_special ≈ particular.u_special
        @test solution.lap.u_lap ≈ lap.u_lap
    end
end
