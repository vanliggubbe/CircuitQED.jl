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
                append!(idxs, last_dof .+ (1 : el.n_aux))
                append!(idxs_el, length(nodes(el)) .+ (1 : el.n_aux))
                last_dof += el.n_aux
            end
            for i in 1 : 3
                K[i][idxs, idxs] .+= K_el[i][idxs_el, idxs_el]
            end
        else
            # nonlinear elements
            # currently only non-biased junctions and SNAILs are allowed
            p_el, q_el, θ_el = nl_response(el, f₀)
            n = length(θ_el)
            p_nl = vcat(p_nl, zeros(n_dof, n))
            q_nl = vcat(p_nl, zeros(n_dof, n))
            p_nl[idxs, end - n + 1 : end] .= p_el[idxs_el, :]
            q_nl[idxs, end - n + 1 : end] .= q_el[idxs_el, :]
            append!(θ_nl, θ_el)
        end
        if el isa Port
            push!(port_Y, inv(unitless(uref(f₀), el.impedance)))
            push!(port_i = first(idxs))
        end
    end
    return ClassicalEOM(f₀, n_dof, K_lin, p_nl, q_nl, θ_nl, port_Y, port_i)
end

struct RFSolver{T <: Real} <: AbstractClassicalSolver
    n_dof :: Int

    RF_lin :: Matrix{Complex{T}}
    RF_inp :: Matrix{Complex{T}}
    RF_p_nl :: Matrix{Complex{T}}

    DC_lin :: Matrix{T}
    DC_p_nl :: Matrix{T}

    q_nl :: Matrix{T}
    θ_nl :: Vector{T}
end

Base.eltype(:: RFSolver{T}) where {T} = T
ndof(solver :: RFSolver) = solver.n_dof * 4 

function rhs!(du, u, p :: Tuple{RFSolver, Function}, t)
    solver, v_i = p
    
    φ_dc = @view u[1 : solver.n_dof] 
    v_dc = @view u[solver.n_dof .+ (1 : solver.n_dof)]
    φ_rf = reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof + (1 : 2 * solver.n_dof)])

    dφ_dc = @view du[1 : solver.n_dof] 
    dv_dc = @view du[solver.n_dof .+ (1 : solver.n_dof)]
    dφ_rf = reinterpret(Complex{eltype(du)}, @view du[2 * solver.n_dof + (1 : 2 * solver.n_dof)])

    # linear dc part
    dφ_dc .= v_dc
    mul!(dv_dc, solver.DC_lin, @view u[1 : (2 * solver.n_dof)])

    # linear rf part
    mul!(dφ_rf, RF_inp, v_i(t))
    mul!(dφ_rf, RF_lin, φ_rf)

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
        axpy!(besselj0(2 * abs_A) * sinb, p_dc, dv_dc)
        # rf part
        axpy!(phs_A * besselj1(2 * abs_A) * cosb, p_rf, dφ_rf)
    end
    return du
end

function jac!(jac, u, p :: Tuple{RFSolver, Function}, t)
    solver, _ = p
    fill!(jac, zero(eltype(jac)))
    n = solver.n_dof

    # linear dc part
    jac[1 : n, n .+ (1 : n)] .= I(n_dof)
    jac[n .+ (1 : n), 1 : 2 * n] .= solver.DC_lin

    # linear rf part
    jac[2 * n .+ (1 : 2 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)] .= reinterpret(Complex{eltype(solver)}, solver.RF_lin)[1 : 2 : end, :]
    jac[2 * n .+ (2 : 2 : 2 * n), 2 * n .+ (2 : 2 : 2 * n)] .= reinterpret(Complex{eltype(solver)}, solver.RF_lin)[1 : 2 : end, :]
    jac[2 * n .+ (2 : 2 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)] .= reinterpret(Complex{eltype(solver)}, solver.RF_lin)[2 : 2 : end, :]
    jac[2 * n .+ (1 : 2 : 2 * n), 2 * n .+ (2 : 2 : 2 * n)] .-= reinterpret(Complex{eltype(solver)}, solver.RF_lin)[2 : 2 : end, :]

    # nonlinear part
    for (i, (p_rf, p_dc, q, θ)) in enumerate(zip(eachcol(solver.RF_p_nl), eachcol(solver.DC_p_nl), eachcol(solver.q_nl), solver.θ_nl))
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
        j0 = besselj0(2 * abs_A)
        j1 = besselj1(2 * abs_A)
        j2 = iszero(abs_A) ? one(eltype(abs_A)) : (j1 / abs_A)

        # dc part
        mul!(view(jac[n .+ (1 : n), 1 : n]), p_dc, q', cosb * j0, one(eltype(jac)))
        mul!(view(jac[n .+ (1 : n), (2 * n) + 1 : 2 : n]), p_dc, q', -sinb * j1 * 2.0 * real(phs_A), one(eltype(jac)))
        mul!(view(jac[n .+ (1 : n), (2 * n) + 2 : 2 : n]), p_dc, q', -sinb * j1 * 2.0 * imag(phs_A), one(eltype(jac)))

        # rf part
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac[2 * n .+ (1 : 2 * n), 1 : n])),
            p_rf, q', -j1 * phs_A * sinb, one(Complex{eltype(jac)})
        )
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac[2 * n .+ (1 : 2 * n), 2 * n .+ (1 : 2 : 2 * n)])),
            p_rf, q', cosb * phs_A * (2.0 * real(phs_A) * j0 - phs_A * j2), one(Complex{eltype(jac)})
        )
        mul!(
            reinterpret(Complex{eltype(jac)}, view(jac[2 * n .+ (1 : 2 * n), 2 * n .+ (2 : 2 : 2 * n)])),
            p_rf, q', im * cosb * phs_A * (-2.0im * imag(phs_A) * j0 + phs_A * j2), one(Complex{eltype(jac)})
        )
        #axpy!(phs_A * besselj1(2 * abs_A) * cosb, p_rf, dφ_rf)
    end

end
