struct ClassicalEOM{T <: Real, F <: Frequency}
    f₀ :: F
    n_dof :: Int
    K_lin :: Tuple{Matrix{T}, Matrix{T}, Matrix{T}}

    M_lin :: Matrix{T}
    M_inp :: Matrix{T}
    M_p_nl :: Matrix{T}

    p_nl :: Matrix{T}
    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
    # current is -K₂ ẍ - K₁ ẋ + K₀ x - p_nl * sin(q_nl' * x - θ_nl)

    input :: SparseMatrixCSC{T, Int}
    #port_Y :: Vector{T}
    #port_i :: Vector{Int}
end

ClassicalEOM(circuit :: Circuit, f₀ :: Frequency) = ClassicalEOM(Float64, circuit, f₀)

function ClassicalEOM(:: Type{T}, circuit :: Circuit, f₀ :: Frequency) where {T <: Real}
    n_dof = ndof(circuit)
    if !isfinite(n_dof)
        error("Cannot deal with infinite number of degrees of freedom")
    end
    K = (
        zeros(T, n_dof, n_dof),
        zeros(T, n_dof, n_dof),
        zeros(T, n_dof, n_dof)
    )
    last_dof = length(circuit.nds) - 1
    p_nl = Matrix{T}(undef, n_dof, 0)
    q_nl = Matrix{T}(undef, n_dof, 0)
    θ_nl = Vector{T}(undef, 0)
    #port_Y = T[]
    #port_i = Int[]
    input = spzeros(T, n_dof, sum(x -> x isa Port, circuit.els))
    i_port = 1
    for el in circuit.els
        idxs = [node_index(circuit, node) for node in nodes(el) if node != ground()]
        idxs_el = [i for (i, node) in enumerate(nodes(el)) if node != ground()]

        if el isa LinearElement
            K_el = response(el, f₀)
            if size(K_el[1], 1) != length(nodes(el))
                # there are auxiliary dofs
                append!(idxs, last_dof .+ (1 : el.n_aux[]))
                append!(idxs_el, length(nodes(el)) .+ (1 : el.n_aux[]))
                last_dof += el.n_aux[]
            end
            for i in 1 : 3
                K[i][idxs, idxs] .+= K_el[i][idxs_el, idxs_el]
            end
        else
            # nonlinear elements
            # currently only non-biased junctions and SNAILs are allowed
            p_el, q_el, θ_el = nl_response(el, f₀)
            n = length(θ_el)
            p_nl = hcat(p_nl, zeros(n_dof, n))
            q_nl = hcat(q_nl, zeros(n_dof, n))
            p_nl[idxs, end - n + 1 : end] .= p_el[idxs_el, :]
            q_nl[idxs, end - n + 1 : end] .= q_el[idxs_el, :]
            append!(θ_nl, θ_el)
        end
        if el isa Port
            #push!(port_Y, inv(unitless(uref(f₀), el.impedance[])))
            #push!(port_i, first(idxs))
            input[first(idxs), i_port] = inv(unitless(uref(f₀), el.impedance[]))
            i_port += 1
        end
    end
    return ClassicalEOM(f₀, n_dof, K, K[3] \ hcat(K[1], -K[2]), -2.0 * K[3] \ input, -K[3] \ p_nl, p_nl, q_nl, θ_nl, input)
end

Base.eltype(:: ClassicalEOM{T}) where {T} = T
ndof(eom :: ClassicalEOM) = 2 * eom.n_dof

function _rhs0!(f, x :: AbstractVector, eom :: ClassicalEOM)
    mul!(f, eom.K_lin[1], x)
    for (p, q, θ) in zip(eachcol(eom.p_nl), eachcol(eom.q_nl), eom.θ_nl)
        axpy!(sin(θ - q' * x), p, f)
    end
    return f
end

function _jac0!(jac, x :: AbstractVector, eom :: ClassicalEOM)
    jac .= eom.K_lin[1]
    for (p, q, θ) in zip(eachcol(eom.p_nl), eachcol(eom.q_nl), eom.θ_nl)
        mul!(jac, p, q', -cos(q' * x - θ), one(eltype(jac)))
    end
    return jac
end


function rhs!(du, u, par :: Tuple{ClassicalEOM, Function}, t)
    (eom, v_i) = par
    φ = @view u[1 : eom.n_dof] 
    v = @view u[eom.n_dof .+ (1 : eom.n_dof)]
    
    dφ = @view du[1 : eom.n_dof] 
    dv = @view du[eom.n_dof .+ (1 : eom.n_dof)]

    # linear part of the EOM
    dφ .= v
    mul!(dv, eom.M_lin, u)
    mul!(dv, eom.M_inp, v_i(t), one(eltype(eom)), one(eltype(dv)))

    for (p, q, θ) in zip(eachcol(eom.M_p_nl), eachcol(eom.q_nl), eom.θ_nl)
        axpy!(sin(q' * φ - θ), p, dv)
    end
    return du
end

function jac!(jac, u, par :: Tuple{ClassicalEOM, Function}, _)
    eom, _ = par
    fill!(jac, zero(eltype(jac)))
    n = eom.n_dof
    φ = @view u[1 : eom.n_dof] 

    # linear part
    @simd for i in 1 : n
        jac[i, n + i] = one(eltype(jac))
    end
    jac[n .+ (1 : n), :] .= eom.M_lin
    for (p, q, θ) in zip(eachcol(eom.M_p_nl), eachcol(eom.q_nl), eom.θ_nl)
        mul!(view(jac, n .+ (1 : n), 1 : n), p, q', cos(q' * φ - θ), one(eltype(jac)))
    end
    return jac
end

"""
    steady_state(eom :: ClassicalEOM; n_newton :: Integer = 10)

Finds stationary solution of classical equations of motion `eom` using Newton method. Optional parameter `n_newton` defines number of iterations.
"""
function steady_state(eom :: ClassicalEOM{T}; n_newton :: Integer = 10) where {T}
    @argcheck n_newton > 1
    x = zeros(T, eom.n_dof)
    y = similar(x)
    f = similar(x)
    jac = zeros(T, eom.n_dof, eom.n_dof)
    for _ in 1 : n_newton
        _rhs0!(f, x, eom)
        _jac0!(jac, x, eom)
        fac = lu!(jac)
        ldiv!(y, fac, f)
        x .-= y
    end
    return [x; zero(x)]
end

"""
    scattering_matrix(eom :: ClassicalEOM, fs; n_newton :: Integer = 10)

For each frequency in `fs`, returns scattering matrix between all ports, assuming the circuit to be in the stationary state.
"""
function scattering_matrix(eom :: ClassicalEOM{T}, fs; n_newton :: Integer = 10) where {T}
    # find stationary state
    x = steady_state(eom; n_newton)
    jac = Matrix{T}(undef, length(x) ÷ 2, length(x) ÷ 2)
    # linearize potential energy around it
    _jac0!(jac, (@view x[1 : length(x) ÷ 2]), eom)

    # construct input and output matrices
    return [
        let ω = unitless(uref(eom.f₀), 2π * f), K0 = jac, K1 = eom.K_lin[2], K2 = eom.K_lin[3], input = -2.0 * eom.input, output = sign.(eom.input), Y = sqrt.(eom.input.nzval)
            Diagonal(Y) * (-1.0I + (-im * ω) * (output' * ((K0 + im * ω * K1 + ω ^ 2 * K2) \ input))) * Diagonal(inv.(Y))
        end
        for f in fs
    ]
end

ODEFunction(eom :: ClassicalEOM) = ODEFunction(rhs!; jac = jac!, jac_prototype = zeros(eltype(eom), ndof(eom), ndof(eom)))

function ODEProblem(eom :: ClassicalEOM, v_i :: Function, tspan :: Tuple{Time, Time}, init = nothing; n_newton :: Integer = 10)
    fun = ODEFunction(eom)
    x₀ = (init isa Nothing ? steady_state(eom; n_newton) : init)
    return ODEProblem(
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
