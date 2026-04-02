struct ClassicalEOM{T <: Real, F <: Frequency}
    f₀ :: F
    n_dof :: Int
    K_lin :: Tuple{Matrix{T}, Matrix{T}, Matrix{T}}
    # current is -K₂ ẍ - K₁ ẋ + K₀ x

    p_nl :: Matrix{T}
    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
    # current is p_nl * sin(q_nl' * x - θ_nl)

    port_Y :: Vector{T}
    port_i :: Vector{Int}
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
    port_Y = T[]
    port_i = Int[]
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
            push!(port_Y, inv(unitless(uref(f₀), el.impedance[])))
            push!(port_i, first(idxs))
        end
    end
    return ClassicalEOM(f₀, n_dof, K, p_nl, q_nl, θ_nl, port_Y, port_i)
end

struct RFSolver{T <: Real} <: AbstractClassicalSolver
    eom :: ClassicalEOM{T}
    n_dof :: Int

    RF_lin :: Matrix{Complex{T}}
    RF_inp :: Matrix{Complex{T}}
    RF_p_nl :: Matrix{Complex{T}}

    DC_lin :: Matrix{T}
    DC_p_nl :: Matrix{T}

    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
end

function RFSolver(eom :: ClassicalEOM, F :: Frequency)
    ω = unitless(uref(eom.f₀), 2π * F)
    input = zeros(eom.n_dof, length(eom.port_i))
    for (j, (Y, i)) in enumerate(zip(eom.port_Y, eom.port_i))
        input[i, j] = 2.0 * Y
    end
    RF = factorize(2.0 * ω * eom.K_lin[3] + im * eom.K_lin[2])
    DC = factorize(eom.K_lin[3])
    return RFSolver(
        eom, eom.n_dof,
        im * (RF \ (ω ^ 2 * eom.K_lin[3] + im * ω * eom.K_lin[2] + eom.K_lin[1])), 
        im * (RF \ input),
        -im * (RF \ eom.p_nl),
        DC \ hcat(eom.K_lin[1], -eom.K_lin[2]),
        -(DC \ eom.p_nl),
        eom.q_nl, eom.θ_nl
    )
end

RFSolver(:: Type{T}, circuit :: Circuit, F :: Frequency, f₀ :: Frequency) where {T <: Real} = RFSolver(ClassicalEOM(T, circuit, f₀), F)
RFSolver(circuit :: Circuit, F :: Frequency, f₀ :: Frequency) = RFSolver(ClassicalEOM(circuit, f₀), F)

Base.eltype(:: RFSolver{T}) where {T} = T
ndof(solver :: RFSolver) = solver.n_dof * 4 

function rhs!(du, u, p :: Tuple{RFSolver, Function}, t)
    solver, v_i = p
    
    φ_dc = @view u[1 : solver.n_dof] 
    v_dc = @view u[solver.n_dof .+ (1 : solver.n_dof)]
    φ_rf = reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof .+ (1 : 2 * solver.n_dof)])

    dφ_dc = @view du[1 : solver.n_dof] 
    dv_dc = @view du[solver.n_dof .+ (1 : solver.n_dof)]
    dφ_rf = reinterpret(Complex{eltype(du)}, @view du[2 * solver.n_dof .+ (1 : 2 * solver.n_dof)])

    # linear dc part
    dφ_dc .= v_dc
    mul!(dv_dc, solver.DC_lin, @view u[1 : (2 * solver.n_dof)])

    # linear rf part
    mul!(dφ_rf, solver.RF_inp, v_i(t))
    mul!(dφ_rf, solver.RF_lin, φ_rf, one(eltype(solver.RF_lin)), one(eltype(dφ_rf)))

    # nonlinear part
    for (p_rf, p_dc, q, θ) in zip(eachcol(solver.RF_p_nl), eachcol(solver.DC_p_nl), eachcol(solver.q_nl), solver.θ_nl)
        A = q' * φ_rf
        b = q' * φ_dc - θ
        #   exp[ i |A| exp(i φ - i Ω t) + i |A| exp(i Ω t - i φ) + i b] / (2 i)
        # - exp[-i |A| exp(i φ - i Ω t) - i |A| exp(i Ω t - i φ) - i b] / (2 i)
        # zero harmonic
        # J₀(2 * |A|) exp(i b) / (2 i) - J₀(2 * |A|) exp(-ib) / (2 i) = J₀
        # first harmonic
        # i exp(i φ) J₁(2 |A|) exp(i b) / (2 i)
        # - i exp(i φ) J₁(-2 |A|) exp(-i b) / (2 i) = 
        # exp(i φ) J₁(2 |A|) [exp(i b) + exp(-i b)] / (2 im)
        
        # precalc
        abs_A = abs(A)
        phs_A = iszero(abs_A) ? one(eltype(A)) : (A / abs_A)
        sinb, cosb = sincos(b)

        # dc part
        axpy!(besselj0(2.0 * abs_A) * sinb, p_dc, dv_dc)
        # rf part
        axpy!(phs_A * besselj1(2.0 * abs_A) * cosb, p_rf, dφ_rf)
    end
    return du
end

function jac!(jac, u, p :: Tuple{RFSolver, Function}, _)
    solver, _ = p
    fill!(jac, zero(eltype(jac)))
    n = solver.n_dof

    φ_dc = @view u[1 : n] 
    φ_rf = reinterpret(Complex{eltype(u)}, @view u[2 * n .+ (1 : 2 * n)])

    # linear dc part
    @simd for i in 1 : n
        jac[i, n + i] = one(eltype(jac))
    end
    jac[n .+ (1 : n), 1 : (2 * n)] .= solver.DC_lin

    # linear rf part
    jac[2 * n .+ (1 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)] .= reinterpret(eltype(solver), solver.RF_lin)
    jac[2 * n .+ (2 : 2 : 2 * n), 2 * n .+ (2 : 2 : 2 * n)] .= @view jac[2 * n .+ (1 : 2 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)]
    jac[2 * n .+ (1 : 2 : 2 * n), 2 * n .+ (2 : 2 : 2 * n)] .-= @view jac[2 * n .+ (2 : 2 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)]

    # nonlinear part
    for (p_rf, p_dc, q, θ) in zip(eachcol(solver.RF_p_nl), eachcol(solver.DC_p_nl), eachcol(solver.q_nl), solver.θ_nl)
        A = q' * φ_rf
        b = q' * φ_dc - θ

        # precalc
        abs_A = abs(A)
        phs_A = iszero(abs_A) ? one(eltype(A)) : (A / abs_A)
        sinb, cosb = sincos(b)
        j0 = besselj0(2 * abs_A)
        j1 = besselj1(2 * abs_A)
        j2 = iszero(abs_A) ? one(eltype(abs_A)) : (j1 / abs_A)

        # dc part
        mul!(view(jac, n .+ (1 : n), 1 : n), p_dc, q', cosb * j0, one(eltype(jac)))
        mul!(view(jac, n .+ (1 : n), (2 * n) .+ (1 : 2 : 2 * n)), p_dc, q', -sinb * j1 * 2.0 * real(phs_A), one(eltype(jac)))
        mul!(view(jac, n .+ (1 : n), (2 * n) .+ (2 : 2 : 2 * n)), p_dc, q', -sinb * j1 * 2.0 * imag(phs_A), one(eltype(jac)))

        # rf part
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac, 2 * n .+ (1 : 2 * n), 1 : n)),
            p_rf, q', -j1 * phs_A * sinb, one(Complex{eltype(jac)})
        )
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac, 2 * n .+ (1 : 2 * n), 2 * n .+ (1 : 2 : 2 * n))),
            p_rf, q', cosb * phs_A * (2.0 * real(phs_A) * j0 - phs_A * j2), one(Complex{eltype(jac)})
        )
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac, 2 * n .+ (1 : 2 * n), 2 * n .+ (2 : 2 : 2 * n))),
            p_rf, q', im * cosb * phs_A * (-2.0im * imag(phs_A) * j0 + phs_A * j2), one(Complex{eltype(jac)})
        )
    end
    return jac
end

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

"""
    steady_state(eom :: ClassicalEOM; n_newton :: Integer = 10)

Finds stationary solution of classical equations of motion `eom` using Newton method. Optional parameter `n_newton` control number of iterations.
"""
function steady_state(eom :: ClassicalEOM{T}; n_newton :: Integer = 10) where {T}
    if n_newton <= 0
        error("Number of Newton iterations must be positive")
    end
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
    return x
end

"""
    steady_state(eom :: ClassicalEOM; n_newton :: Integer = 10)

Finds stationary solution of classical equations of motion `rf` using Newton method.
"""
steady_state(rf :: RFSolver{T}; n_newton :: Integer = 10) where {T} = [steady_state(rf.eom; n_newton); zeros(T, 3 * rf.n_dof)]

"""
    scattering_matrix(eom :: ClassicalEOM, fs; n_newton :: Integer = 10)

For each frequency in `fs`, returns scattering matrix between all ports, assuming the circuit to be in the stationary state.
"""
function scattering_matrix(eom :: ClassicalEOM{T}, fs; n_newton :: Integer = 10) where {T}
    # find stationary state
    x = steady_state(eom; n_newton)
    jac = Matrix{T}(undef, length(x), length(x))
    # linearize potential energy around it
    _jac0!(jac, x, eom)

    # construct input and output matrices
    input = spzeros(T, eom.n_dof, length(eom.port_i))
    output = spzeros(T, eom.n_dof, length(eom.port_i))
    for (j, (Y, i)) in enumerate(zip(eom.port_Y, eom.port_i))
        input[i, j] = -2.0 * Y
        output[i, j] = one(eltype(T))
    end
    return [
        let ω = unitless(uref(eom.f₀), 2π * f), K0 = jac, K1 = eom.K_lin[2], K2 = eom.K_lin[3], input = input, output = output
            -1.0I + (-im * ω) * (output' * ((K0 + im * ω * K1 + ω ^ 2 * K2) \ input))
        end
        for f in fs
    ]
end

scattering_matrix(rf :: RFSolver, fs; n_newton :: Integer = 10) = scattering_matrix(rf.eom, fs; n_newton)

ODEFunction(rf :: RFSolver) = ODEFunction(rhs!; jac = jac!, jac_prototype = zeros(eltype(rf), ndof(rf), ndof(rf)))

function ODEProblem(rf :: RFSolver, v_i :: Function, tspan :: Tuple{Time, Time}, init = nothing; n_newton :: Integer = 10)
    fun = ODEFunction(rf)
    x₀ = (init isa Nothing ? steady_state(rf; n_newton) : init)
    return ODEProblem(
        fun, x₀, (
            unitless(uref(rf.eom.f₀), tspan[1]),
            unitless(uref(rf.eom.f₀), tspan[2])
        ), (
            rf, let v_i = v_i, f₀ = rf.eom.f₀;
                t -> [unitless(uref(f₀), x) for x in v_i(unitof(Time, uref(f₀)) * t)]
            end
        )
    )
end
