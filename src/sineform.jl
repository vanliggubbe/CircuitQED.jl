struct AffineSineForm{T, TA <: AbstractMatrix{T}, TB <: AbstractMatrix{T}, TP <: AbstractMatrix{T}, TQ <: AbstractMatrix{T}, TT <: AbstractVector{T}}
    # y = A * x + B * u + p_nl * sin.(q_nl' * x .- θ_nl)
    A :: TA
    B :: TB
    p_nl :: TP
    q_nl :: TQ
    θ_nl :: TT
end

function _eval!(y, form :: AffineSineForm, x)
    mul!(y, form.A, x)
    for (p, q, θ) in zip(eachcol(form.p_nl), eachcol(form.q_nl), form.θ_nl)
        axpy!(sin(q' * x - θ), p, y)
    end
    return y
end

function _eval!(y, form :: AffineSineForm, x, u)
    _eval!(y, form, x)
    mul!(y, form.B, u, one(eltype(y)), one(eltype(y)))
    return y
end

function _jac!(jac, form :: AffineSineForm, x)
    jac .= form.A
    for (p, q, θ) in zip(eachcol(form.p_nl), eachcol(form.q_nl), form.θ_nl)
        mul!(jac, p, q', cos(q' * x - θ), one(eltype(jac)))
    end
    return jac
end

function (f :: AffineSineForm{T})(x) where {T}
    y = Vector{promote_type(eltype(x), T)}(undef, size(f.A, 1))
    return _eval!(y, f, x)
end

function (f :: AffineSineForm{T})(x, u) where {T}
    y = Vector{promote_type(eltype(x), eltype(u), T)}(undef, size(f.A, 1))
    return _eval!(y, f, x, u)
end


