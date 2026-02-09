import ProgressLogging
using Plots

include("../src/mod.jl")

Δ_scan(p, t) = p.Δ₀ + p.Δ_rate * t
 p = (;
    S=7.5,
    Δ₀=-4.0,
    Δ_rate=1.5,
    Δ = Δ_scan, 
    τ_R=2e-4,
)
dτ = 0.05

u0 = lle_hss(p.S, p.Δ₀, 1024; noise=1e-4)

t = 0.0:0.05:50
LLE_raman = SemilinearPDE(
    StandartLLE(;
         N̂=LLERamanNonlinearPart() # Raman nonlinearity
    ),
    u0,
    t[end],
    dτ,
    p
);

u2_raman = solve_rkip(LLE_raman; saveat=t, abstol=1e-2, reltol=1e-4);
heatmap(abs2.(u2_raman))






