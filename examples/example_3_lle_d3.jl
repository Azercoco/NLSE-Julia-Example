import ProgressLogging
include("../src/mod.jl")
using Plots

Δ_scan(p, t) = p.Δ₀ + p.Δ_rate * t
p = (;
    S=4.5,
    Δ₀=-2.0,
    Δ_rate=0.4,
    Δ=Δ_scan,
    d₃=0.08,
)
dτ = 0.025

u0 = lle_hss(p.S, p.Δ₀, 2^12; noise=1e-2)

t = 0.0:0.25:50
LLE_d3 = SemilinearPDE(
    StandartLLE(;
        D̂=(ω, p) -> -1 - 1im * ω^2 + 1im * ω^3 * p.d₃
    ),
    u0,
    t[end],
    dτ,
    p
);

u2_d3 = solve_rkip(LLE_d3; saveat=t, abstol=1e-4, reltol=1e-8);
heatmap(abs2.(u2_d3))

begin
    # How to use DSP to compute a PSD 
    using FFTW
    using DSP

    prgd = welch_pgram(u2_d3[:, end], 2^9; fs=inv(dτ), window=bartlett_hann)

    plot(
        prgd.freq |> fftshift, #fftshift ensures the freqences are sorted
        10 .* log10.(prgd.power) |> fftshift,
        xlims=(-10, 10)
    )
end

