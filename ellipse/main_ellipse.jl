using Gridap
using GridapGmsh
using LinearAlgebra
using GridapMakie, CairoMakie
using Arpack
using IterativeSolvers
using Printf

model = GmshDiscreteModel(joinpath(@__DIR__, "geo-ellipse.msh"))

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")

g(x) = 0.0
Ug = TrialFESpace(V0, g)

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

a1(u, v) = ∫(∇(u)⋅∇(v))*dΩ
a2(u, v) = ∫(u*v)*dΩ
K = assemble_matrix(a1, Ug, V0)
M = assemble_matrix(a2, Ug, V0)

λ, ϕ = eigs(K, M; nev=3, which=:SM)

k = sqrt(λ[1])
F(x) = -k*cos(x[1]+x[2])

eig_funcs = [FEFunction(Ug, ϕ[:, i]) for i in 1:3]

# fig = Figure(size = (1650, 480), fontsize = 18)
# for i in 1:3
#     ax = Axis(fig[1, 2i-1],
#               title = @sprintf("Mode %d\nλ_%d = %.4f", i, i, λ[i]),
#               xlabel = "x", ylabel = i == 1 ? "y" : "",
#               aspect = DataAspect())
#     plt = plot!(ax, Ω, eig_funcs[i], colormap = :RdBu)
#     Colorbar(fig[1, 2i], plt, width = 12)
# end
# save(joinpath(@__DIR__, "eigenfunctions.png"), fig) 

c_1 = sum(∫(F*eig_funcs[1])*dΩ)

A = K - k^2*M
l(v) = ∫(F*v)dΩ - c_1*∫(eig_funcs[1]*v)dΩ
b = assemble_vector(l, V0)

u_special_vec = zeros(length(b))
u_special_vec, history = gmres!(
    u_special_vec, A, b; restart=50, maxiter=25000, reltol=1e-6, log=true
)
history

u_special = FEFunction(Ug, u_special_vec)

# fig, ax, plt = plot(Ω, u_special)
# Colorbar(fig[1, 2], plt)
# fig

a_1 = -c_1 / (2*k^2)
u_final_vec = u_special_vec+a_1*ϕ[:, 1]
u_final = FEFunction(Ug, u_final_vec)

fig, ax, plt = plot(Ω, u_final)
Colorbar(fig[1, 2], plt)
fig
