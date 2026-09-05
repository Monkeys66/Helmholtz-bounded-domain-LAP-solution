# Assemble the finite element matrices.
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


# Assemble the source and its k-derivative.
function assemble_source_loads(fem, source, source_k, k)
    dΩ = fem.dΩ

    F(x) = source(x, k)
    F_k(x) = source_k(x, k)

    l(v) = ∫(F * v) * dΩ
    l_k(v) = ∫(F_k * v) * dΩ

    b = assemble_vector(l, fem.V0)
    b_der = assemble_vector(l_k, fem.V0)

    return (; b, b_der)
end
