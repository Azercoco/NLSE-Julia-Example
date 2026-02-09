import ProgressLogging
include("../src/mod.jl")
using Plots

Δ_scan(p, t) = p.Δ₀ + p.Δ_rate * t
p = (;
    S=3.5,
    Δ₀=-2.0,
    Δ_rate=0.2,
    Δ = Δ_scan, 
)
dτ = 0.05

u0 = lle_hss(p.S, p.Δ₀, 1024; noise=1e-2)

t = 0.0:0.25:100
LLE_scan = SemilinearPDE(
    StandartLLE(),
    u0,
    t[end],
    dτ,
    p
);

u2_scan = solve_rkip(LLE_scan; saveat=t, abstol=1e-4, reltol=1e-8);
heatmap(abs2.(u2_scan))

