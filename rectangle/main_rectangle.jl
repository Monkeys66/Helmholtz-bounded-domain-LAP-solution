using Gridap
using GridapGmsh
using LinearAlgebra
using GridapMakie, CairoMakie
using Arpack
using IterativeSolvers

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

k = sqrt(5.0)
F(x) = -x[1]*x[2]

eig_func1 = FEFunction(Ug, ϕ[:, 1])
eig_func2 = FEFunction(Ug, ϕ[:, 2])
eig_func3 = FEFunction(Ug, ϕ[:, 3])

# fig, ax, plt = plot(Ω, eig_func3)
# Colorbar(fig[1,2], plt)
# fig

c_1 = sum(∫(F*eig_func2)*dΩ) 
c_2 = sum(∫(F*eig_func3)*dΩ) 

A = K - k^2*M
l(v) = ∫(F*v)dΩ - c_1*∫(eig_func2*v)dΩ - c_2*∫(eig_func3*v)dΩ
b = assemble_vector(l, V0)

u_special_vec = zeros(length(b))
u_special_vec, history = gmres!(
    u_special_vec, A, b; restart=50, maxiter=25000, reltol=1e-10, log=true
)
history

u_special = FEFunction(Ug, u_special_vec)

fig, ax, plt = plot(Ω, u_special)
Colorbar(fig[1, 2], plt)
fig

a_1 = c_1/10
a_2 = c_2/10

u_final_vec = u_special_vec + a_1*ϕ[:,2] + a_2*ϕ[:,3]
u_final = FEFunction(Ug,u_final_vec)

fig, ax, plt = plot(Ω, u_final)
Colorbar(fig[1,2], plt)
fig



