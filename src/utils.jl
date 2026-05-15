const Mayhaps = Union{Nothing, T} where {T}

@inline uref(f :: Frequency) = (2π * f, 1.0u"ħ", 2.0u"q", 1.0u"k")

function _rowspace_basis(A :: AbstractMatrix{T}) where {T <: Real}
    if size(A, 2) == 0
        return zeros(T, 0, 0)
    end
    decomp = svd(Matrix(A))
    σmax = isempty(decomp.S) ? zero(T) : maximum(decomp.S)
    tol = max(size(A)...) * eps(T) * σmax
    r = count(>(tol), decomp.S)
    return Matrix(decomp.V[:, 1 : r])
end

function _range_null_basis(A :: AbstractMatrix{T}) where {T <: Real}
    if size(A, 1) == 0
        return zeros(T, 0, 0), zeros(T, 0, 0)
    end
    decomp = svd(Matrix(A))
    σmax = isempty(decomp.S) ? zero(T) : maximum(decomp.S)
    tol = max(size(A)...) * eps(T) * σmax
    r = count(>(tol), decomp.S)
    V = Matrix(decomp.V)
    return V[:, 1 : r], V[:, r + 1 : end]
end


function inf_toeplitz_fun(left :: Int, coeff :: Vector{T}, fun :: Function; atol :: Real = zero(real(T)), rtol :: Real = iszero(atol) ? sqrt(eps(real(T))) : zero(real(T))) where {T <: Number}
    @argcheck atol >= zero(eltype(atol))
    @argcheck rtol >= zero(eltype(rtol))

    right = left + length(coeff) - 1
    # Work in a floating complex type
    R = float(real(T))
    C = Complex{R}

    # Exact zero case: exp(0) = I
     
    if all(iszero, coeff)
        return T <: Real ? (0, 0, [one(R)]) : (0, 0, [one(C)])
    end

    N = 4
    while N < max(abs(left) + 1, abs(right) + 1)
        N *= 2
    end

    r = zeros(C, N)
    r_prev = C[]
    while true
        fill!(r, zero(C))
        idxs = (left : right) .+ ((left : right) .< 0) * length(r) .+ 1
        r[idxs] .= coeff
        fft!(r)
        @inbounds @simd for i in eachindex(r)
            r[i] = fun(r[i])
        end
        ifft!(r)
        if !isempty(r_prev)
            k = length(r_prev) ÷ 2
            tmp = sum(abs, r_prev[1 : k] - r[1 : k]) + sum(abs, r_prev[k .+ (1 : k)] - r[3 * k .+ (1 : k)]) + sum(abs, r[k .+ (1 : 2 * k)])
            if tmp < sum(abs, r) * rtol
                break
            end
        end
        r_prev = copy(r)
        append!(r, similar(r))
    end
    circshift!(r, length(r) ÷ 2)
    return (-length(r) ÷ 2, T <: Real ? real(r) : r)
end
