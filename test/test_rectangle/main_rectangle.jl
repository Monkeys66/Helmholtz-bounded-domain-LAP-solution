using Gridap
using GridapGmsh
using LinearAlgebra
using GridapMakie, CairoMakie
using Arpack
using IterativeSolvers
using Printf

model = GmshDiscreteModel(joinpath(@__DIR__, "geo-rectangle.msh"))

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

λ, ϕ = eigs(K, M; nev=4, which=:SM)

k = sqrt(λ[2])
F(x) = -x[1]*x[2]

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

c_1 = sum(∫(F*eig_funcs[2])*dΩ)
c_2 = sum(∫(F*eig_funcs[3])*dΩ)

A = K - k^2*M
l(v) = ∫(F*v)dΩ - c_1*∫(eig_funcs[2]*v)dΩ - c_2*∫(eig_funcs[3]*v)dΩ
b = assemble_vector(l, V0)

u_special_vec = zeros(length(b))
u_special_vec, history = gmres!(
    u_special_vec, A, b; restart=50, maxiter=25000, reltol=1e-10, log=true
)
history

u_special = FEFunction(Ug, u_special_vec)

a_1 = c_1/10
a_2 = c_2/10

u_final_vec = u_special_vec + a_1*ϕ[:, 2] + a_2*ϕ[:, 3]
u_final = FEFunction(Ug, u_final_vec)

# fig = Figure(; size=(820, 440), fontsize=16)
# ax = Axis(
#     fig[1, 1];
#     title="special solution-rectangle case",
#     xlabel="x",
#     ylabel="y",
#     aspect=DataAspect(),
# )
# plt = plot!(ax, Ω, u_special; colormap=:RdBu)
# Colorbar(fig[1, 2], plt)
# save(joinpath(@__DIR__, "u_special_rectangle.png"), fig; px_per_unit=3)

# fig = Figure(; size=(820, 440), fontsize=16)
# ax = Axis(
#     fig[1, 1];
#     title="LAP solution-rectangle case",
#     xlabel="x",
#     ylabel="y",
#     aspect=DataAspect(),
# )
# plt = plot!(ax, Ω, u_final; colormap=:RdBu)
# Colorbar(fig[1, 2], plt)
# save(joinpath(@__DIR__, "u_final_rectangle.png"), fig; px_per_unit=3)

# fig

