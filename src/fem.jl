#assemble the finite element matrices 
function assemble_fem(meshfile; order=1, degree=2, boundary_tag="boundary")
    model = GmshDiscreteModel(meshfile)

    reffe = ReferenceFE(lagrangian, Float64, order)

    V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=boundary_tag)

    g(x) = 0.0
    U = TrialFESpace(V0, g)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    a_K(u, v) = ∫(∇(u) ⋅ ∇(v)) * dΩ
    a_M(u, v) = ∫(u * v) * dΩ

    K = assemble_matrix(a_K, U, V0)
    M = assemble_matrix(a_M, U, V0)

    return (; model, V0, U, Ω, dΩ, K, M)
end

#compute the eigenvalus and eigenfuctions
function compute_eigenpairs(K, M; nev, which=:SM)
    λ, Φ = eigs(K, M; nev=nev, which=which)

    return (; λ, Φ)
end

#print the eigenvalues
function show_eigens(λ)
    len = length(λ)
    for i in 1:len
        @printf("Eingenvalue λ_%d = %.4f\n", i, i, λ[i])
    end
end


#compute the particular solution 
function compute_parSolution(K, M, dΩ, V0, λ, Φ, F, start_index, end_index)
    k = sqrt(λ[start_index])
    len = end_index - start_index + 1

    eigs = λ[start_index:end_index]
    eig_funcs = [FEFunction(V0, Φ[:, i]) for i in start_index:end_index]

    coefficients_kernel = [sum(∫(F * eig_funcs[i]) * dΩ) for i in 1:len]
    sum_coefficients = sum(coefficients_kernel)

    for i in 1:len
        if !isapprox(coefficients_kernel[i], 0.0; atol=1e-4)
            @warn "F is not orthogonal to the kernel space!"
        end
    end

    A = K - k^2 * M
    l(v) = ∫(F * v) * dΩ - sum(∫(coefficients_kernel[i] * eig_funcs[i] * v) * dΩ for i in 1:len)
    b = assemble_vector(l, V0)

    u_special_vec = zeros(length(b))

    u_special_vec, history =
        gmres!(u_special_vec, A, b; restart=50, maxiter=25000, reltol=1e-10, log=true)

    return u_special_vec
end


