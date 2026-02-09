using LinearAlgebra
using SciMLOperators: AbstractSciMLOperator
using UnPack


"""
    DiagonalOperator{T}

Spectral linear operator represented by a diagonal matrix.

This operator is designed for Fourier-space representations of linear
PDE terms, where differentiation becomes multiplication by a diagonal
symbol.

It is compatible with `SplitODEProblem` and SciML operator interfaces.
"""
struct DiagonalOperator{T,uType<:AbstractVector{T}} <: AbstractSciMLOperator{T}
    diag::uType
end

"""
    exp(D::DiagonalOperator, t)

Compute the matrix exponential `exp(t * D)`.

Since `D` is diagonal, the exponential is applied elementwise to the
diagonal entries.

Returns a new `DiagonalOperator`.
"""
function LinearAlgebra.exp(D::DiagonalOperator{T,uType}, t) where {T,uType}
    exp_kernel = copy(D.diag)
    exp_kernel .*= t
    exp_kernel .= exp.(exp_kernel)
    return DiagonalOperator{T,uType}(exp_kernel)
end

"""
    (D::DiagonalOperator)(du, u, _, _, _)

Apply the diagonal operator in-place.

This matches the calling convention expected by SciML split ODE solvers.
"""
(D::DiagonalOperator)(du, u, _, _, _) = du .= D.diag .* u


"""
Abstract supertype for nonlinear PDE terms.

A nonlinear part must implement a call overload with signature:

`
nl(du, u, p, t)
`

Optionally, it may provide:
- `get_cache(nl::AbstractNonlinearPart, u, p)` for preallocated work arrays. Returns `nothing` by default.
- time- or parameter-dependent coefficients via `get_var(nl::AbstractNonlinearPart, var, p, t)`.

`p` will contains the equations parameters as well as
- `p.τ`: array spatial variable.
- `p.freq`: array, corresponding spectral variable.
- `p.bifft_plan`: `BiFFTPlan`, for preallocated FFTs and IFFTs.


Cache initialized with `get_cache(nl, u, p)` will be passed through p 
and acessible with `p.cache`.
"""
abstract type AbstractNonlinearPart end


"""
    get_cache(nl, u, p)

Return a cache object for the nonlinear term.

By default, no cache is used.
"""
get_cache(::AbstractNonlinearPart, _, _) = nothing

"""
    get_var(nl, var, p, t)

Resolve a parameter that may be either a constant or a function of `(p, t)`.

This allows nonlinear terms to use time-dependent coefficients
without branching at every call site.
"""
get_var(::AbstractNonlinearPart, var::Real, _, _) = var
get_var(::AbstractNonlinearPart, var::Function, p, t) = var(p, t)

"""
    NonlinearPartWrapper

Wraps a nonlinear operator that acts in physical space, while the solver
state is stored in spectral space.

This wrapper:
1. Transforms `u` from spectral to physical space,
2. Applies the nonlinear operator in physical space,
3. Transforms the result back to spectral space.

This pattern is used to reduce the number of FFTs by half in the RKIP implementation.
"""
struct NonlinearPartWrapper{F<:AbstractNonlinearPart,tmpType,cacheType}
    f::F
    tmp::tmpType
    cache::cacheType
end

"""
Apply the wrapped nonlinear operator.

Input and output are in spectral space.
Internally, the computation is performed in physical space using FFTs.
"""
@fastmath @inline function (N̂::NonlinearPartWrapper)(du, u, p, t)
    bifft_plan = p.bifft_plan
    fft_tmp = bifft_plan.fft_tmp

    @unpack tmp, f, cache = N̂

    copyto!(fft_tmp, u)
    bifft_plan.ifft_plan * fft_tmp

    copyto!(tmp, fft_tmp)

    if isnothing(cache)
        f(du, tmp, p, t)
    else
        f(du, tmp, (; cache, p...), t)
    end

    copyto!(fft_tmp, du)
    bifft_plan.fft_plan * fft_tmp
    copyto!(du, fft_tmp)

    return du
end