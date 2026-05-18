struct AffineSineForm{T <: Real}
    # y = A * x + B * u + p_nl * sin.(q_nl' * x .- θ_nl)
    A :: Matrix{T}
    B :: Matrix{T}
    p_nl :: Matrix{T}
    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
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

function (f :: AffineSineForm)(x, u)
    y = Vector{promote_type(eltype(x), eltype(u))}(undef, size(f.A, 1))
    return _eval!(y, f, x, u)
end

struct ClassicalEOM{T <: Real, F <: Frequency}
    f₀ :: F

    # Minimal explicit state equation:
    # ẋ = dyn.A * x + dyn.B * v_in(t) +
    #     dyn.p_nl * sin.(dyn.q_nl' * x .- dyn.θ_nl).
    dyn :: AffineSineForm{T}

    # Port voltage readout:
    # v_out = out * ẋ - v_in(t)
    out :: Matrix{T}

    # Maps Port element names to rows of the output form.
    port_index :: Dict{Symbol, Int}
end

ClassicalEOM(circuit :: Circuit, f₀ :: Frequency) = ClassicalEOM(Float64, circuit, f₀)

function _add_response!(K, K_el, idxs, idxs_el)
    for (i, i_el) in zip(idxs, idxs_el)
        for (j, j_el) in zip(idxs, idxs_el)
            K[i, j] += K_el[i_el, j_el]
        end
    end
    return K
end

function _add_nl_response!(dst, src, idxs, idxs_el, cols)
    for (i, i_el) in zip(idxs, idxs_el)
        dst[i, cols] .+= src[i_el, :]
    end
    return dst
end

function ClassicalEOM(:: Type{T}, circuit :: Circuit, f₀ :: Frequency) where {T <: Real}
    n_raw = ndof(circuit)
    @argcheck isfinite(n_raw) "Number of degrees of freedom must be finite"

    # dynamic matrices related to inductive, resistive, and capacitive elements
    K = (zeros(T, n_raw, n_raw), zeros(T, n_raw, n_raw), zeros(T, n_raw, n_raw))
    # nonlinear terms related to the Josephson junction
    p_nl = Matrix{T}(undef, n_raw, 0)
    q_nl = Matrix{T}(undef, n_raw, 0)
    θ_nl = Vector{T}(undef, 0)
    # ports
    input = zeros(T, n_raw, sum(x -> x isa Port, circuit.els))
    port_voltage = zeros(T, size(input, 2), n_raw)
    port_index = Dict{Symbol, Int}()

    i_port = 1
    for el in circuit.els
        idxs, idxs_el = _linear_idxs(circuit, el)
        if el isa LinearElement
            if el isa Short
                continue
            end
            K_el = response(el, f₀)
            @argcheck all(size(K_i, 1) == length(coordinates(el)) && size(K_i, 2) == length(coordinates(el)) for K_i in K_el)
            for i in 1 : 3
                _add_response!(K[i], K_el[i], idxs, idxs_el)
            end
        else
            p_el, q_el, θ_el = nl_response(el, f₀)
            n = length(θ_el)
            p_nl = hcat(p_nl, zeros(T, n_raw, n))
            q_nl = hcat(q_nl, zeros(T, n_raw, n))
            cols = (size(p_nl, 2) - n + 1) : size(p_nl, 2)
            _add_nl_response!(p_nl, p_el, idxs, idxs_el, cols)
            _add_nl_response!(q_nl, q_el, idxs, idxs_el, cols)
            append!(θ_nl, θ_el)
        end
        if el isa Port
            if !isempty(idxs)
                input[first(idxs), i_port] = inv(unitless(uref(f₀), el.impedance[]))
                port_voltage[i_port, first(idxs)] = one(T)
            end
            port_index[el.name] = i_port
            i_port += 1
        end
    end

    # find bases of kernels and ranges
    n_φ = n_raw # number of flux DOF
    n_v = n_raw # number of voltage DOF
    
    local dyn
    local out
    if rank(K[3]) < size(K[3], 1)
        # we have damped modes
        # or even algebraic conditions
        U, S, V = svd(K[3])
        rk = rank(Diagonal(S))

        # transform to coordinates which diagonalize massive part
        K3 = Diagonal(S)
        K2 = U' * K[2] * V
        K1 = U' * K[1] * V
        q_nl = V' * q_nl
        p_nl = U' * p_nl
        input = U' * input
        port_voltage = port_voltage * V

        if rank(K2[rk + 1 : end, rk + 1 : end]) == n_φ - rk
            # only damped modes
            n_v = rk
            dyn = AffineSineForm(
                zeros(T, n_φ + rk, n_φ + rk),
                zeros(T, n_φ + rk, size(input, 2)),
                zeros(T, n_φ + rk, size(p_nl, 2)),
                [q_nl; zeros(T, rk, size(q_nl, 2))],
                θ_nl
            )

            # deal with damped modes
            M2 = factorize(K2[rk + 1 : end, rk + 1 : end])
            dyn.A[1 : rk, n_φ + 1 : n_φ + rk] = I(rk)
            dyn.A[rk + 1 : n_φ, 1 : n_φ] = (M2 \ K1[rk + 1 : n_φ, :])
            dyn.A[rk + 1 : n_φ, n_φ + 1 : n_φ + rk] = -(M2 \ K2[rk + 1 : end, 1 : rk])
            dyn.B[rk + 1 : n_φ, :] = 2 * (M2 \ input[rk + 1 : end, :])
            dyn.p_nl[rk + 1 : n_φ, :] = -(M2 \ p_nl[rk + 1 : n_φ, :])

            # deal with massive modes
            M3 = K3[1 : rk, 1 : rk] \ K2[1 : rk, :]
            dyn.A[n_φ + 1 : end, 1 : n_φ] = (K3[1 : rk, 1 : rk] \ K1[1 : rk, :])
            dyn.A[n_φ + 1 : end, :] -= M3 * dyn.A[1 : n_φ, :]
            dyn.B[n_φ + 1 : end, :] = 2 * (K3[1 : rk, 1 : rk] \ input[1 : rk, :])
            dyn.B[n_φ + 1 : end, :] -= M3 * dyn.B[1 : n_φ, :]
            dyn.p_nl[n_φ + 1 : end, :] = -K3[1 : rk, 1 : rk] \ p_nl[1 : rk, :]
            dyn.p_nl[n_φ + 1 : end, :] -= M3 * dyn.p_nl[1 : n_φ, :]
            out = hcat(port_voltage, zeros(T, size(input, 2), rk))
        else
            # there are algebraic constraints
            error("Circuits with algebraic constraints are not supported yet")
        end
    else
        dyn = AffineSineForm(
            [
                zeros(T, n_φ, n_φ)  I(n_φ);
                (K[3] \ K[1])       (-(K[3] \ K[2]))
            ],
            [
                zeros(T, n_φ, size(input, 2));
                2 * (K[3] \ input)
            ],
            [zeros(T, n_φ, size(p_nl, 2)); -K[3] \ p_nl],
            [q_nl; zeros(T, n_v, size(q_nl, 2))],
            θ_nl
        )
        out = hcat(port_voltage, zeros(T, size(input, 2), n_v))
    end

    cokernel, kernel = domain_basis([dyn.A; dyn.q_nl'])
    if size(kernel, 1) > 0
        # there are cyclic degrees of freedom
        dyn = AffineSineForm(
            cokernel' * dyn.A * cokernel,
            cokernel' * dyn.B,
            cokernel' * dyn.p_nl,
            cokernel' * dyn.q_nl,
            dyn.θ_nl
        )
        out = out * cokernel
    end
    return ClassicalEOM(f₀, dyn, out, port_index)
end

Base.eltype(:: ClassicalEOM{T}) where {T} = T
ndof(eom :: ClassicalEOM) = size(eom.dyn.A, 1)

_rhs0!(f, x :: AbstractVector, eom :: ClassicalEOM) = _eval!(f, eom.dyn, x)
_jac0!(jac, x :: AbstractVector, eom :: ClassicalEOM) = _jac!(jac, eom.dyn, x)

function rhs!(du, u, par :: Tuple{ClassicalEOM, Function}, t)
    eom, v_i = par
    return _eval!(du, eom.dyn, u, v_i(t))
end

function jac!(jac, u, par :: Tuple{ClassicalEOM, Function}, _)
    eom, _ = par
    return _jac!(jac, eom.dyn, u)
end

output_voltage(eom :: ClassicalEOM, u, v_in) = (eom.out * eom.dyn(u, v_in) - v_in)

"""
    steady_state(eom :: ClassicalEOM; n_newton :: Integer = 10)

Finds a stationary solution of the explicit state equation using Newton method. Optional parameter `n_newton` defines number of iterations.
"""
function steady_state(eom :: ClassicalEOM{T}; n_newton :: Integer = 10) where {T}
    @argcheck n_newton > 1
    x = zeros(T, ndof(eom))
    y = similar(x)
    f = similar(x)
    jac = zeros(T, ndof(eom), ndof(eom))
    for _ in 1 : n_newton
        _rhs0!(f, x, eom)
        _jac0!(jac, x, eom)
        fac = lu!(jac)
        ldiv!(y, fac, f)
        x .-= y
    end
    return x
end

function scattering_matrix(eom :: ClassicalEOM, fs; n_newton :: Integer = 10)
    steady = steady_state(eom; n_newton)
    jac = schur(_jac0!(zeros(eltype, ndof(eom), ndof(eom)), steady, eom))
    ref = uref(eom.f₀)
    S = Matrix{Complex{eltype(eom)}}[]
    for f in fs
        ω = 2π * unitless(ref, f)
        push!(S, (-im * ω) * eom.out * jac.Z * ((-im * ω - jac.T) \ (jac.Z' * eom.dyn.B)) - I)
    end
    return S
end


ODEFunction(eom :: ClassicalEOM) = ODEFunction(rhs!; jac = jac!, jac_prototype = zeros(eltype(eom), ndof(eom), ndof(eom)))

function ODEProblem(eom :: ClassicalEOM, v_i :: Function, tspan :: Tuple{Time, Time}, init = nothing; n_newton :: Integer = 10)
    fun = ODEFunction(eom)
    x₀ = (init isa Nothing ? steady_state(eom; n_newton) : init)
    return ODEProblem{true}(
        fun, x₀, (
            unitless(uref(eom.f₀), tspan[1]),
            unitless(uref(eom.f₀), tspan[2])
        ), (
            eom, let v_i = v_i, f₀ = eom.f₀;
                t -> [unitless(uref(f₀), x) for x in v_i(unitof(Time, uref(f₀)) * t)]
            end
        )
    )
end
