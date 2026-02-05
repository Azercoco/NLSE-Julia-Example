import ProgressLogging
using Plots

include("../src/mod.jl")

# FFTW.set_provider!("mkl")
FFTW.set_num_threads(14)

Δ_scan(p, t) = p.Δ₀ + p.Δ_rate * t
 
p = (;
    S=5.5,
    Δ₀=-4.0,
    Δ_rate=0.15,
    τ_R=2e-4,
    Δ = Δ_scan, 
)

u0 = lle_hss(p.S, p.Δ₀, 1024; noise=1e-4)
dτ = 0.05
t = 0.0:2.5:300

LLE_scan = SemilinearPDE(
    StandartLLE(;
         N̂=LLERamanNonlinearPart() # Raman nonlinearity
    ),
    u0,
    t[end],
    dτ,
    p
);

u2_scan = solve_rkip(LLE_scan; saveat=t);

heatmap(abs2.(u2_scan))

# u2_scan = solve_rkip(LLE_scan)





