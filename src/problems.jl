using SciMLBase
using DiffEqBase
using OrdinaryDiffEqRKIP: RKIP
using DiffEqCallbacks


"""
    SeminlinearPDEOperators

Container for the operators defining a semilinear PDE.

# Parameters

- `D̂`: Function defining the linear operator symbol in Fourier space.
        Signature: `(ω, p) -> scalar`
- `N̂::AbstractNonlinearPart`: Nonlinear operator acting in physical space.
"""
struct SeminlinearPDEOperators{F1<:Function,F2<:AbstractNonlinearPart}
    D̂::F1
    N̂::F2
end

"""
    init(op, u0, p)

Initialize the linear and nonlinear operators for a semilinear PDE.

Returns:
- `D̂`: DiagonalOperator for the linear part (spectral space),
- `N̂`: NonlinearPartWrapper for the nonlinear part.
"""
function init(op::SeminlinearPDEOperators, u0, p)
    D̂ = DiagonalOperator(
        convert(typeof(u0), op.D̂.(2π .* p.freq, Ref(p)))
    )

    N̂ = NonlinearPartWrapper(
        op.N̂,
        similar(u0),
        get_cache(op.N̂, u0, p)
    )

    return D̂, N̂
end

"""
    SemilinearPDE

High-level representation of a semilinear PDE solved using a
Fourier pseudospectral method and operator splitting.

This struct bundles together FFT plans, operators, grids, and the
underlying `SplitODEProblem`.


"""
struct SemilinearPDE{tPlan<:BiFFTPlan,tOps<:SeminlinearPDEOperators,tTau,tFreq,tProblem}
    bifft_plan::tPlan
    ops::tOps
    τ::tTau
    freq::tFreq
    problem::tProblem
end


"""
    SemilinearPDE(op, u0, tspan, dτ, p)

Construct a semilinear PDE problem from operators and initial data.

# Parameters

- `u0` :  Complex vector of the initial condition, given *in physical space*.
- `tspan` : Tuple or real containing the tntegration bounds, if only a real is given, assumed to be `(0, tspan)`
- `dτ` : Real, 
- `p` : Named Tuple, storing the PDE parameters


Internally, the problem is transformed to spectral space to reduce the number of FFTs.
The resulting object is ready to be solved with time integrators
such as `RKIP`.

Returns a `SemilinearPDE` instance.
"""
function SemilinearPDE(op::SeminlinearPDEOperators,
    u0::AbstractVector{Complex{T}},
    tspan,
    dτ::Real,
    p
) where {T<:Real}
    # @assert (length(u0) % number_of_variables == 0) "Numbers of DOF = $(length(u)) is not an integer multiple of the number of variables = $(number_of_variables)"

    nb_pts = length(u0) #÷ number_of_variables
    domain_size = dτ * nb_pts

    bifft_plan = BiFFTPlan(u0)

    τ = similar(u0, T)
    freq = similar(u0, T)

    p_ext = (; τ, bifft_plan, freq, p...)

    τ .= range(-domain_size / 2.0, domain_size / 2.0, nb_pts)
    freq .= FFTW.fftfreq(nb_pts, inv(dτ))

    D̂, N̂ = init(op, u0, p_ext)

    tspan = isa(tspan, Real) ? (0, tspan) : tspan #if no initial value is given, use zero

    problem = SplitODEProblem(
        SplitFunction{true}(D̂, N̂),
        fft(bifft_plan, u0),
        tspan,
        p_ext,
    )

    return SemilinearPDE{typeof(bifft_plan),typeof(op),typeof(τ),typeof(freq),typeof(problem)}(bifft_plan, op, τ, freq, problem)
end
"""
    solve_rkip!(res, pde; kwargs...)

Solve a `SemilinearPDE` using the RKIP integrator.

The solution is returned in **physical space**.

If `saveat` is provided, `res` must be a matrix whose columns correspond
to saved time snapshots.

Other parameter will be passed to `OrdinaryDiffEq.solve`
"""
function solve_rkip!(res,
    pde::SemilinearPDE;
    alg=RKIP(1e-4, 1.0),
    reltol=1e-6,
    abstol=1e-10,
    progress=true,
    saveat=nothing,
    progress_steps=100,
    save_everystep=false,
    dense=false,
    kwargs...
)

    if isnothing(saveat)
        sol = DiffEqBase.solve(pde.problem, alg; reltol, abstol, save_everystep, progress, progress_steps, dense, kwargs...)
        ifft!(pde.bifft_plan, res, sol.u[end])
    else
        @assert size(res)[2] == length(saveat) "The ouput with size $(size(res)) is not compatible with the size of saveat $(length(saveat))"
        index = Ref(1) #Reference to be mutable outside the loop

        # in-place writing to avoid using too much memory
        save_at_index(integrator) =
            let
                ifft!(pde.bifft_plan, view(res, :, index[]), integrator.u)
                index[] += 1
            end

        sol = DiffEqBase.solve(pde.problem, alg; reltol, abstol, progress, save_everystep, progress_steps, dense, callback=PresetTimeCallback(
                saveat,
                save_at_index;
                save_positions=(false, false),
            ), kwargs...)
    end
    return res
end


"""
    solve_rkip(pde; kwargs...)

Convenience wrapper allocating the output array automatically.

Returns the solution in physical space.
"""
solve_rkip(pde::SemilinearPDE; saveat=nothing, kwargs...) = solve_rkip!(
    isnothing(saveat) ? similar(pde.problem.u0) : similar(pde.problem.u0,
        (length(pde.problem.u0), length(saveat))
    ), pde; saveat, kwargs...
)
