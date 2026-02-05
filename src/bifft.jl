using FFTW
"""
    BiFFTPlan{Vc, FFT_type, IFFT_type}

Bidirectional FFT plan used to efficiently transform vectors between
physical space and spectral (Fourier) space **in-place**.

This structure owns:
- a temporary buffer `fft_tmp` used for FFT/IFFT operations,
- a pre-planned forward FFT (`fft_plan`),
- a pre-planned inverse FFT (`ifft_plan`).

The plan is created once and reused to avoid repeated allocations and
planning costs, which is critical in time-stepping PDE solvers.

# Fields
- `fft_tmp`: Temporary array used internally for FFT operations.
- `fft_plan`: FFTW plan for forward FFT.
- `ifft_plan`: FFTW plan for inverse FFT.
"""
struct BiFFTPlan{Vc,FFT_type,IFFT_type}
    fft_tmp::Vc
    fft_plan::FFT_type
    ifft_plan::IFFT_type
end

"""
    BiFFTPlan(u)

Construct a `BiFFTPlan` compatible with the array `u`.

If `u` is real-valued, the internal buffer is promoted to a complex type,
since FFTs operate in the complex domain.

The returned plan can be reused for forward and inverse FFTs on arrays
with the same size and element type as `u`.
"""
function BiFFTPlan(u)
    T = eltype(u)
    T = (T <: Real) ? Complex{T} : T

    tmp = similar(u, T)

    return BiFFTPlan(
        tmp,
        plan_fft!(tmp),
        plan_ifft!(tmp)
    )

end


"""
    fft_op!(plan, res, u, fft_op)

Internal helper applying a pre-planned FFT or IFFT operator.

Copies `u` into the plan's temporary buffer, applies the FFT operator,
and writes the result into `res`.

This function is allocation-free assuming `res` is preallocated.
"""
function fft_op!(plan::BiFFTPlan, res, u, fft_op)
    copyto!(plan.fft_tmp, u)
    fft_op * plan.fft_tmp
    copyto!(res, plan.fft_tmp)
    return res
end

"""
Compute the forward in-place FFT of `u` into `res` using `plan`.
"""
fft!(plan::BiFFTPlan, res, u) = fft_op!(plan, res, u, plan.fft_plan)

"""
Compute the inverse in-place FFT of `u` into `res` using `plan`.
"""
ifft!(plan::BiFFTPlan, res, u) = fft_op!(plan, res, u, plan.ifft_plan)

"""
Compute the forward in-place FFT of `u` using `plan`, overwriting `u`'
"""
fft!(plan::BiFFTPlan, u) = fft!(plan, u, u)
"""
Compute the inverse in-place FFT of `u` using `plan`, overwriting `u`'
"""
ifft!(plan::BiFFTPlan, u) = ifft!(plan, u, u)


"""
Out-of-place FFT `u`.

Allocates a new array for the result.
"""
fft(plan::BiFFTPlan, u) = fft!(plan, similar(u), u)

"""
Out-of-place IFFT `u`.

Allocates a new array for the result.
"""
ifft(plan::BiFFTPlan, u) = ifft!(plan, similar(u), u)