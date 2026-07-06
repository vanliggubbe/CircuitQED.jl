struct ClassicalEOM{T <: Real, U <: Tuple{Vararg{Quantity}}, FD <: AffineSineForm{T}, FO <: AffineSineForm{T}, FJ <: AffineSineForm{T}}
    units :: U
    circuit :: Circuit

    # Minimal explicit state equation:
    # ẋ = dyn.A * x + dyn.B * v_in(t) +
    #     dyn.p_nl * sin.(dyn.q_nl' * x .- dyn.θ_nl).
    dyn :: FD

    # Port voltage readout:
    # v_out = out.A * x + out.B * v_in(t) +
    #         out.p_nl * sin.(out.q_nl' * x .- dyn.θ_nl).
    out :: FO

    # Voltage accross each of the resistors
    # used to calculate Joule losses
    jls :: FJ
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
    input = zeros(T, n_raw, length(circuit.port_list))
    port_voltage = zeros(T, size(input, 2), n_raw)
    resistor_voltage = zeros(T, length(circuit.resistor_index), n_raw)
    units = uref(f₀)

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
            i_port = circuit.port_index[el.name]
            if !isempty(idxs)
                input[first(idxs), i_port] = inv(unitless(units, el.impedance[]))
                port_voltage[i_port, first(idxs)] = one(T)
            end
        elseif el isa Resistor
            # TODO make it better
            i_res = circuit.resistor_index[el.name]
            r = sqrt(unitless(units, el.resistance[]))
            resistor_voltage[i_res, first(idxs)] = -inv(r)
            resistor_voltage[i_res,  last(idxs)] =  inv(r)
        end
    end

    # find bases of kernels and ranges
    n_φ = n_raw # number of flux DOF
    n_v = n_raw # number of voltage DOF
    
    local dyn
    local out
    local jls
    rk = rank(K[3])
    if rk < size(K[3], 1)
        # we have damped modes
        # or even algebraic conditions
        U, S, V = svd(K[3])

        # transform to coordinates which diagonalize massive part
        K3 = Diagonal(S)[1 : rk, 1 : rk]
        K2 = U' * K[2] * V
        K1 = U' * K[1] * V
        q_nl = V' * q_nl
        p_nl = U' * p_nl
        input = U' * input
        port_voltage        = port_voltage * V
        resistor_voltage    = resistor_voltage * V

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
            dyn.B[rk + 1 : n_φ, :] = 2 * (M2 \ input[rk + 1 : n_φ, :])
            dyn.p_nl[rk + 1 : n_φ, :] = -(M2 \ p_nl[rk + 1 : n_φ, :])

            # deal with massive modes
            dyn.A[n_φ + 1 : end, n_φ + 1 : end] = -(K3 \ K2[1 : rk, 1 : rk])
            dyn.A[n_φ + 1 : end, 1 : n_φ] = (K3 \ K1[1 : rk, :])
            dyn.B[n_φ + 1 : end, :] = 2 * (K3 \ input[1 : rk, :])
            dyn.p_nl[n_φ + 1 : end, :] = -(K3 \ p_nl[1 : rk, :])
            # whatever comes from damped modes
            dyn.A[n_φ + 1 : end, :] -= (K3 \ K2[1 : rk, rk + 1 : n_φ]) * dyn.A[rk + 1 : n_φ, :]
            dyn.B[n_φ + 1 : end, :] -= (K3 \ K2[1 : rk, rk + 1 : n_φ]) * dyn.B[rk + 1 : n_φ, :]
            dyn.p_nl[n_φ + 1 : end, :] -= (K3 \ K2[1 : rk, rk + 1 : n_φ]) * dyn.p_nl[rk + 1 : n_φ, :]

            tmp = hcat(port_voltage, zeros(T, size(input, 2), rk))
            out = AffineSineForm(tmp * dyn.A, tmp * dyn.B - I, tmp * dyn.p_nl, dyn.q_nl, dyn.θ_nl)

            tmp = hcat(resistor_voltage, zeros(T, size(resistor_voltage, 1), rk))
            jls = AffineSineForm(tmp * dyn.A, tmp * dyn.B, tmp * dyn.p_nl, dyn.q_nl, dyn.θ_nl)
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
        tmp = hcat(port_voltage, zeros(T, size(input, 2), n_v))
        out = AffineSineForm(tmp * dyn.A, tmp * dyn.B - I, tmp * dyn.p_nl, dyn.q_nl, dyn.θ_nl)

        tmp = hcat(resistor_voltage, zeros(T, size(resistor_voltage, 1), n_v))
        jls = AffineSineForm(tmp * dyn.A, tmp * dyn.B, tmp * dyn.p_nl, dyn.q_nl, dyn.θ_nl)
    end

    cokernel, kernel = domain_basis([dyn.A; dyn.q_nl'])
    if size(kernel, 2) > 0
        # there are cyclic degrees of freedom
        dyn = AffineSineForm(
            cokernel' * dyn.A * cokernel,
            cokernel' * dyn.B,
            cokernel' * dyn.p_nl,
            cokernel' * dyn.q_nl,
            dyn.θ_nl
        )
        out = AffineSineForm(
            out.A * cokernel,
            out.B,
            out.p_nl,
            dyn.q_nl,
            dyn.θ_nl
        )
        jls = AffineSineForm(
            jls.A * cokernel,
            jls.B,
            jls.p_nl,
            dyn.q_nl,
            dyn.θ_nl
        )
    end
    return ClassicalEOM(units, circuit, dyn, out, jls)
end

Base.eltype(:: ClassicalEOM{T}) where {T} = T
ndof(eom :: ClassicalEOM) = size(eom.dyn.A, 1)
_port_admittances(eom :: ClassicalEOM) = (unitless(eom.units, Y) for Y in _port_admittances(eom.circuit))

_rhs0!(f, x :: AbstractVector, eom :: ClassicalEOM) = _eval!(f, eom.dyn, x)
_jac0!(jac, x :: AbstractVector, eom :: ClassicalEOM) = _jac!(jac, eom.dyn, x)

_unitless_input_vector(units, v_i, port_index :: AbstractDict{Symbol, Int}, dim :: Int) = sparsevec(
    Dict(
        let i = port_index[p], u = u, units = units
            @argcheck u isa Voltage
            i => unitless(units, u)
        end
        for (p, u) in v_i
    ), dim
)

function rhs!(du, u, par :: Tuple{ClassicalEOM, Function}, t)
    eom, v_i = par
    return _eval!(du, eom.dyn, u, v_i(t))
end

function jac!(jac, u, par :: Tuple{ClassicalEOM, Function}, _)
    eom, _ = par
    return _jac!(jac, eom.dyn, u)
end


"""
    steady_state(eom :: ClassicalEOM, v_i = (); init = zeros(eltype(eom), ndof(eom)), reltol = sqrt(eps(eltype(eom))), kwargs...)

Find a stationary solution for constant input voltages. `v_i` may be any iterable of `port => voltage` pairs.
Extra keyword arguments are passed to `NonlinearSolve.solve`. Throws an error if the final relative residual norm exceeds `reltol`.
"""
function steady_state(eom :: ClassicalEOM{T}, v_i = (p.name => 0.0u"V" for p in eom.circuit.port_list); init = zeros(T, ndof(eom)), reltol :: Real = sqrt(eps(T)), kwargs...) where {T}
    v_in = _unitless_input_vector(eom.units, v_i, eom.circuit.port_index, length(eom.circuit.port_list))

    function f!(res, x, _)
        return _eval!(res, eom.dyn, x, v_in)
    end
    function j!(jac, x, _)
        return _jac!(jac, eom.dyn, x)
    end

    fun = NonlinearFunction(f!; jac = j!)
    prob = NonlinearSolve.NonlinearProblem(fun, init, nothing)
    sol = NonlinearSolve.solve(prob; kwargs...)
    x = sol.u
    res = similar(x)
    _eval!(res, eom.dyn, x, v_in)
    relres = norm(res) / max(norm(x), one(T))
    if !(relres < reltol)
        error("Steady state solve did not converge: relative residual $(relres) exceeds $(reltol)")
    end
    return x
end

"""
    scattering_matrix(eom :: ClassicalEOM, fs; kwargs...)

Returns a vector of scattering matrices for each frequency in `fs`, assuming the circuit is in the steady state. Keyword arguments are passed to `steady_state`.
"""
function scattering_matrix(eom :: ClassicalEOM, fs; kwargs...)
    steady = steady_state(eom; kwargs...)
    fac = schur(_jac0!(zeros(eltype(eom), ndof(eom), ndof(eom)), steady, eom))
    ports = port_names(eom.circuit)
    out = _jac!(zeros(eltype(eom), length(ports), ndof(eom)), eom.out, steady) * fac.Z
    right = fac.Z' * eom.dyn.B
    return [
        let ω = 2π * unitless(eom.units, f), ports = ports, L = fac.T, left = out, right = right, B = eom.out.B
            NamedArray(left * ((-im * ω * I - L) \ right) + B; names = (ports, ports))
        end for f in fs
    ]
end

"""
    scattering_matrix(eom :: ClassicalEOM, f :: Unitful.Frequency; kwargs...)

Returns a scattering matrix for frequency `f`, assuming the circuit is in the steady state. Keyword arguments are passed to `steady_state`.
"""
function scattering_matrix(eom :: ClassicalEOM, f :: Frequency; kwargs...)
    steady = steady_state(eom; kwargs...)
    fac = schur(_jac0!(zeros(eltype(eom), ndof(eom), ndof(eom)), steady, eom))
    ports = port_names(eom.circuit)
    out = _jac!(zeros(eltype(eom), length(ports), ndof(eom)), eom.out, steady) * fac.Z
    right = fac.Z' * eom.dyn.B
    return let ω = 2π * unitless(eom.units, f), ports = ports, L = fac.T, left = out, right = right, B = eom.out.B
        NamedArray(left * ((-im * ω * I - L) \ right) + B; names = (ports, ports))
    end
end


ODEFunction(eom :: ClassicalEOM) = ODEFunction(rhs!; jac = jac!, jac_prototype = zeros(eltype(eom), ndof(eom), ndof(eom)))

unitless_input(
    units :: Tuple{Vararg{Quantity}}, 
    v_i :: Function, 
    circuit :: Circuit
) = let v_i = v_i, port_index = circuit.port_index, t₀ = unitof(Time, units), units = units, dim = length(circuit.port_list)
    t -> _unitless_input_vector(units, v_i(t * t₀), port_index, dim)
end

function ODEProblem(eom :: ClassicalEOM, v_i :: Function, tspan :: Tuple{Time, Time}, init = nothing)
    fun = ODEFunction(eom)
    x₀ = (init isa Nothing ? steady_state(eom) : init)
    return ODEProblem{true}(
        fun, x₀, (
            unitless(eom.units, tspan[1]),
            unitless(eom.units, tspan[2])
        ), (
            eom, unitless_input(eom.units, v_i, eom.circuit)
        )
    )
end

_output_voltage(eom :: ClassicalEOM, u, v_in) = eom.out(u, v_in)
_output_current(eom :: ClassicalEOM, u, v_in) = collect(_port_admittances(eom)) .* (_output_voltage(eom, u, v_in) .+ 2 .* v_in)
_joule_losses(eom :: ClassicalEOM, u, v_in) = eom.jls(u, v_in) .^ 2

"""
    output_voltage(sol :: ODESolution, t :: Time)

Return output voltages at all ports at time `t`.
"""
function output_voltage(sol :: ODESolution, t :: Time)
    (eom, v_i) = sol.prob.p
    t_ul = unitless(eom.units, t)
    V₀ = unitof(u"V", eom.units)
    return NamedArray(_output_voltage(eom, sol(t_ul), v_i(t_ul)) * V₀; names = (port_names(eom.circuit), ))
end

_port_name(x) = Symbol(x)
_port_name(x :: Symbol) = x
_port_name(x :: Port) = x.name

"""
    output_voltage(sol :: ODESolution, port, t :: Time)

Return output voltage at `port` at time `t`.
"""
function output_voltage(sol :: ODESolution, port, t :: Time)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.port_index[_port_name(port)]
    t_ul = unitless(eom.units, t)
    return _output_voltage(eom, sol(t_ul), v_i(t_ul))[i] * unitof(u"V", eom.units)
end

"""
    output_voltage(sol :: ODESolution, port, ts)

Return output voltage at `port` for each time in `ts`.
"""
function output_voltage(sol :: ODESolution, port, ts)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.port_index[_port_name(port)]
    return [
        let t_ul = unitless(eom.units, t);
            _output_voltage(eom, sol(t_ul), v_i(t_ul))[i]
        end for t in ts
    ] * unitof(u"V", eom.units)
end

"""
    output_current(sol :: ODESolution, t :: Time)

Return output currents at all ports at time `t`.
"""
function output_current(sol :: ODESolution, t :: Time)
    (eom, v_i) = sol.prob.p
    t_ul = unitless(eom.units, t)
    I₀ = unitof(u"A", eom.units)
    return NamedArray(_output_current(eom, sol(t_ul), v_i(t_ul)) * I₀; names = (port_names(eom.circuit), ))
end

"""
    output_current(sol :: ODESolution, port, t :: Time)

Return output current at `port` at time `t`.
"""
function output_current(sol :: ODESolution, port, t :: Time)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.port_index[_port_name(port)]
    t_ul = unitless(eom.units, t)
    return _output_current(eom, sol(t_ul), v_i(t_ul))[i] * unitof(u"A", eom.units)
end

"""
    output_current(sol :: ODESolution, port, ts)

Return output current at `port` for each time in `ts`.
"""
function output_current(sol :: ODESolution, port, ts)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.port_index[_port_name(port)]
    return [
        let t_ul = unitless(eom.units, t);
            _output_current(eom, sol(t_ul), v_i(t_ul))[i]
        end for t in ts
    ] * unitof(u"A", eom.units)
end

"""
    joule_losses(sol :: ODESolution, t :: Time)

Return power of Joule losses at all resistors at time `t`.
"""
function joule_losses(sol :: ODESolution, t :: Time)
    (eom, v_i) = sol.prob.p
    t_ul = unitless(eom.units, t)
    P₀ = unitof(u"W", eom.units)
    return NamedArray(_joule_losses(eom, sol(t_ul), v_i(t_ul)) * P₀; names = (_resistor_names(eom.circuit), ))
end


_resistor_name(x) = Symbol(x)
_resistor_name(x :: Symbol) = x
_resistor_name(x :: Resistor) = x.name

"""
    joule_losses(sol :: ODESolution, resistor, t :: Time)

Return power of Joule losses at `resistor` at time `t`.
"""
function joule_losses(sol :: ODESolution, res, t :: Time)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.resistor_index[_resistor_name(res)]
    t_ul = unitless(eom.units, t)
    return _joule_losses(eom, sol(t_ul), v_i(t_ul))[i] * unitof(u"W", eom.units)
end

"""
    joule_losses(sol :: ODESolution, resistor, ts)

Return power of Joule losses at `resistor` for each time in `ts`.
"""
function joule_losses(sol :: ODESolution, res, ts)
    (eom, v_i) = sol.prob.p
    i = eom.circuit.resistor_index[_resistor_name(res)]
    return [
        let t_ul = unitless(eom.units, t);
            _joule_losses(eom, sol(t_ul), v_i(t_ul))[i]
        end for t in ts
    ] * unitof(u"W", eom.units)
end
