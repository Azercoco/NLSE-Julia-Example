import ProgressLogging
include("../src/mod.jl")
using Plots

p = (; S=1.1, Δ=1.1);

u0 = lle_hss(p.S, p.Δ, 2^12; noise=1e-2)

dτ = 0.05
LLE = SemilinearPDE(
    StandartLLE(),
    u0,
    300.0,
    dτ,
    p
);

u = solve_rkip(LLE; abstol=1e-4, reltol=1e-8)
plot(
    LLE.τ,
    abs2.(u),
)

t = 0.0:5:300
u2 = solve_rkip(LLE; saveat=t, abstol=1e-4, reltol=1e-8);

heatmap(abs2.(u2))