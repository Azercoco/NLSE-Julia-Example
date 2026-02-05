using Polynomials
using LoopVectorization

"""
    lle_hss(S, Δ)

Computes the lower Homegenous Steady State (HSS) field solution of the Lugiato–Lefever equation.

# Parameters

- `S` : Real, normalized driving.
- `Δ` : Real, normalized detuning.

# Optional Parameters

- `size` : Create a uniform Vector of length `size` filled with the HSS solution.
- `T` : Convert the vector to type `Complex{T}`. Default is `Float64`.

# Keyword Parameters
- `noise` : Add a random complex vector with amplitude `noise` and a random uniform phase to the HSS. Default is `0.0`. 
Only works if `size` is specified.


"""
function lle_hss(S, Δ)
    coeffs = [-S^2, (1 + Δ^2), -2Δ, 1]
    p = Polynomial(coeffs)
    r = roots(p)
    Y = minimum(real.(filter(isreal, r)))
    return S / (1 + 1im * (Δ - Y))
end

lle_hss(S, Δ, size, T, ; noise=0.0) = fill(Complex{T}(lle_hss(S, Δ)), size) .+ Complex{T}.(noise .* cis.(2π .* rand(size)))
lle_hss(S, Δ, size; noise) = lle_hss(S, Δ, size, Float64, ; noise)



"""
    StandartLLE(; D̂, N̂)

Construct a `SeminlinearPDEOperators` object corresponding to the
**normalized Lugiato–Lefever equation (LLE)**.

This helper provides sensible defaults for both the linear spectral
operator and the nonlinear physical-space operator, matching the
standard normalized form used in nonlinear optics and Kerr resonator
models.

# Mathematical form

The normalized LLE reads:

∂ₜ u = -(1 + iΔ)u + i|u|²u - i∂²_τ u + S

where:
- `Δ` is the detuning,
- `S` is the pump field,
- `τ` is the fast-time coordinate.

The equation is split into:
- a **linear part** handled in Fourier space,
- a **nonlinear part** handled in physical space.

# Keyword Arguments

- `D̂`: Function defining the linear operator symbol in Fourier space.

  *Signature:* `(ω, p) -> scalar`

  *Default:*
  `
  (ω, p) -> -1 - 1im * ω^2
  ` for anomalous dispersion

-  `N̂::AbstractNonlinearPart`: Nonlinear operator acting in physical space.
  Default: `LLENonlinearPart()`

"""
StandartLLE(;
    D̂=(ω, p) -> -1 - 1im * ω^2,
    N̂=LLENonlinearPart()
) = SeminlinearPDEOperators(
    D̂,
    N̂,
)

"""
    LLENonlinearPart <: AbstractNonlinearPart

Nonlinear operator for the normalized Lugiato–Lefever equation (LLE).

This operator represents the **local Kerr nonlinearity, detuning,
and coherent driving** acting in physical space.

This operator computes: N(u) = `i(|u|² - Δ)u + S`

where:
- `Δ` is the detuning parameter,
- `S` is the pump field amplitude.

Both `Δ` and `S` may be:
- constants, or
- functions of `(p, t)` for time-dependent forcing.

# Required parameters in `p`

The parameter container `p` must define:
- `p.Δ` — detuning (real or function),
- `p.S` — pump field (complex or function).

"""
struct LLENonlinearPart <: AbstractNonlinearPart

end



function @fastmath (lle::LLENonlinearPart)(du, u, p, t, cache)
    Δ = get_var(lle, p.Δ, p, t)
    S = get_var(lle, p.S, p, t)

    @. du = 1im * (abs2.(u) - Δ) * u + S
end




