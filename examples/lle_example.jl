import ProgressLogging
using Plots

include("../src/mod.jl")

##=

p = (; S=1.1, Δ=1.1);

u0 = lle_hss(p.S, p.Δ, 2^10; noise=1e-2)
dτ = 0.05

LLE = SemilinearPDE(
    StandartLLE(),
    u0,
    300.0,
    dτ,
    p
);

##=

u = solve(LLE; progress_steps=1)
plot(
    LLE.τ,
    abs2.(u),
)
##=

saveat = 0.0:5:300
u2 = solve_rkip(LLE; saveat);

heatmap(abs2.(u2))