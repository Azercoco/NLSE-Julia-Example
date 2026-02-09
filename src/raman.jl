
"""
    LLERamanNonlinearPart <: AbstractNonlinearPart

Nonlinear operator for the normalized Lugiato–Lefever equation (LLE).

This operator represents the **local Kerr nonlinearity, detuning,
and coherent driving and linearized Raman nonlinearity** acting in physical space.

This operator computes: N(u) = `i(|u|² - Δ)u + S + (i τ_R * ∂_τ |u|²) u ` 

where:
- `Δ` is the detuning parameter,
- `S` is the pump field amplitude.
- `τ_R` is the normalized Raman response time. 

Both `Δ` and `S` may be:
- constants, or
- functions of `(p, t)` for time-dependent forcing.

# Required parameters in `p`

The parameter container `p` must define:
- `p.Δ` — detuning (real or function),
- `p.S` — pump field (complex or function).
- `p.τ_R` normalized Raman response time (real). 

"""
struct LLERamanNonlinearPart <: AbstractNonlinearPart end

# For Raman, we need to implement a cache to store intermediate computation
function get_cache(::LLERamanNonlinearPart, u0, p)
    abs2_tmp = similar(u0) # pre-allocate an array for intermediate storage
    raman_fourier_resp = convert(
        typeof(u0),
        @. -1im * 2π * p.freq * p.τ_R
    ) # pre-computed the Rama kernel in spectral space
    return (; abs2_tmp, raman_fourier_resp) # cache are supplied as a NamedTuple
end


@fastmath function (lle_raman::LLERamanNonlinearPart)(du, u, p, t)
    S = get_var(lle_raman, p.S, p, t) # fetch eventually time varying function
    Δ = get_var(lle_raman, p.Δ, p, t) # fetch eventually time varying function

    @unpack abs2_tmp, raman_fourier_resp = p.cache # We can acces the cache using p.cache

    abs2_tmp .= abs2.(u) # used form temporary storage

    @. du = 1im * (abs2_tmp - Δ) * u + S # LLE

    fft!(p.bifft_plan, abs2_tmp) 
    @. abs2_tmp = raman_fourier_resp * abs2_tmp
    ifft!(p.bifft_plan, abs2_tmp)

    @. du += 1im * abs2_tmp * u

    return du
end