struct AffineSineForm{T <: Real}
    # y = A*x + B*u + p_nl*sin.(q_nl'*x .- θ_nl)
    A :: Matrix{T}
    B :: Matrix{T}
    p_nl :: Matrix{T}
    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
end

struct ClassicalEOM{T <: Real, F <: Frequency}
    f₀ :: F

    # Minimal explicit state equation:
    # ẋ = dynamics.A*x + dynamics.B*v_in(t) +
    #      dynamics.p_nl*sin.(dynamics.q_nl'*x .- dynamics.θ_nl).
    dynamics :: AffineSineForm{T}

    # Port voltage readout:
    # v_out = output.A*x + output.B*v_in(t) +
    #         output.p_nl*sin.(output.q_nl'*x .- output.θ_nl).
    output :: AffineSineForm{T}

    # Maps Port element names to rows of the output form.
    port_index :: Dict{Symbol, Int}

    # Raw circuit phases are coordinate_basis * reduced_phase_coordinates.
    coordinate_basis :: Matrix{T}
end

ClassicalEOM(circuit :: Circuit, f₀ :: Frequency) = ClassicalEOM(Float64, circuit, f₀)

function _linear_idxs(circuit, el)
    idxs = Int[]
    idxs_el = Int[]
    for (i, node) in enumerate(coordinates(el))
        idx = node_index(circuit, node)
        if idx != 0
            push!(idxs, idx)
            push!(idxs_el, i)
        end
    end
    return idxs, idxs_el
end

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

function _eval!(y, form :: AffineSineForm, x, u)
    mul!(y, form.A, x)
    mul!(y, form.B, u, one(eltype(y)), one(eltype(y)))
    for (p, q, θ) in zip(eachcol(form.p_nl), eachcol(form.q_nl), form.θ_nl)
        axpy!(sin(q' * x - θ), p, y)
    end
    return y
end

function _jac!(jac, form :: AffineSineForm, x)
    jac .= form.A
    for (p, q, θ) in zip(eachcol(form.p_nl), eachcol(form.q_nl), form.θ_nl)
        mul!(jac, p, q', cos(q' * x - θ), one(eltype(jac)))
    end
    return jac
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

    # Remove phase directions that do not affect equations, nonlinear phases,
    # inputs, or observable port voltages.
    R = _rowspace_basis([K[1]; K[2]; K[3]; q_nl'; input'; port_voltage])
    K0 = R' * K[1] * R
    K1 = R' * K[2] * R
    K2 = R' * K[3] * R
    p_red = R' * p_nl
    q_red = R' * q_nl
    input_red = R' * input
    port_voltage_red = port_voltage * R

    U, N = _range_null_basis(K2)
    n_inertial = size(U, 2)
    n_massless = size(N, 2)
    n_coord = size(R, 2)
    n_state = 2 * n_inertial + n_massless
    n_input = size(input, 2)
    n_terms = length(θ_nl)

    # State layout is [a; ȧ; b], where reduced phases are ψ = U*a + N*b.
    phase_from_state = hcat(U, zeros(T, n_coord, n_inertial), N)
    velocity_from_adot = hcat(zeros(T, n_coord, n_inertial), U, zeros(T, n_coord, n_massless))
    q_state = phase_from_state' * q_red

    bA = zeros(T, n_massless, n_state)
    bB = zeros(T, n_massless, n_input)
    bP = zeros(T, n_massless, n_terms)
    if n_massless > 0
        Dnn = N' * K1 * N
        @argcheck size(_rowspace_basis(Dnn), 2) == n_massless "ClassicalEOM has algebraic constraints after gauge reduction; explicit ODE reduction is not implemented for this circuit"
        bA .= Dnn \ (N' * K0 * phase_from_state - N' * K1 * velocity_from_adot)
        bB .= Dnn \ (-2.0 * N' * input_red)
        bP .= Dnn \ (-N' * p_red)
    end

    aA = zeros(T, n_inertial, n_state)
    aB = zeros(T, n_inertial, n_input)
    aP = zeros(T, n_inertial, n_terms)
    if n_inertial > 0
        Mii = U' * K2 * U
        aA .= Mii \ (U' * K0 * phase_from_state - U' * K1 * velocity_from_adot - U' * K1 * N * bA)
        aB .= Mii \ (-2.0 * U' * input_red - U' * K1 * N * bB)
        aP .= Mii \ (-U' * p_red - U' * K1 * N * bP)
    end

    A = zeros(T, n_state, n_state)
    B = zeros(T, n_state, n_input)
    P = zeros(T, n_state, n_terms)
    if n_inertial > 0
        A[1 : n_inertial, n_inertial .+ (1 : n_inertial)] .= Matrix{T}(I, n_inertial, n_inertial)
        A[n_inertial .+ (1 : n_inertial), :] .= aA
        B[n_inertial .+ (1 : n_inertial), :] .= aB
        P[n_inertial .+ (1 : n_inertial), :] .= aP
    end
    if n_massless > 0
        rows = 2 * n_inertial .+ (1 : n_massless)
        A[rows, :] .= bA
        B[rows, :] .= bB
        P[rows, :] .= bP
    end

    # Port voltages are linear combinations of phase velocities. For massless
    # coordinates, phase velocity is itself obtained from the explicit RHS.
    velocity_from_state = velocity_from_adot + N * bA
    C = port_voltage_red * velocity_from_state
    D = port_voltage_red * N * bB
    P_out = port_voltage_red * N * bP

    dynamics = AffineSineForm(A, B, P, q_state, θ_nl)
    output = AffineSineForm(C, D, P_out, q_state, θ_nl)
    return ClassicalEOM(f₀, dynamics, output, port_index, R)
end

Base.eltype(:: ClassicalEOM{T}) where {T} = T
ndof(eom :: ClassicalEOM) = size(eom.dynamics.A, 1)

function _rhs0!(f, x :: AbstractVector, eom :: ClassicalEOM)
    return _eval!(f, eom.dynamics, x, zeros(eltype(eom), size(eom.dynamics.B, 2)))
end

function _jac0!(jac, x :: AbstractVector, eom :: ClassicalEOM)
    return _jac!(jac, eom.dynamics, x)
end

function rhs!(du, u, par :: Tuple{ClassicalEOM, Function}, t)
    eom, v_i = par
    return _eval!(du, eom.dynamics, u, v_i(t))
end

function jac!(jac, u, par :: Tuple{ClassicalEOM, Function}, _)
    eom, _ = par
    return _jac!(jac, eom.dynamics, u)
end

function output_voltage(eom :: ClassicalEOM, u, v_in)
    y = similar(u, size(eom.output.A, 1))
    return _eval!(y, eom.output, u, v_in)
end

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
    error("scattering_matrix for first-order ClassicalEOM is not implemented yet")
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
