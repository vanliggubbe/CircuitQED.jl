struct RFSolver{T <: Real} <: AbstractClassicalSolver
    eom :: ClassicalEOM{T}
    ω :: T
    n_fourier :: Int

    #RF_lin :: Matrix{Complex{T}}
    #RF_inp :: Matrix{Complex{T}}
    #RF_p_nl :: Matrix{Complex{T}}
    RF_lin  :: Array{Complex{T}, 3}
    RF_inp  :: Array{Complex{T}, 3}
    RF_p_nl :: Array{Complex{T}, 3}

    DC_lin  :: Matrix{T}
    DC_p_nl :: Matrix{T}
end

function Base.getproperty(rf :: RFSolver, name :: Symbol)
    if name == :q_nl
        return rf.eom.q_nl
    elseif name == :θ_nl
        return rf.eom.θ_nl
    elseif name == :n_dof
        return rf.eom.n_dof
    else
        return getfield(rf, name)
    end
end

function RFSolver(eom :: ClassicalEOM{T}, F :: Frequency, n_fourier :: Int = 1) where {T}
    @argcheck n_fourier > 0
    @argcheck F > zero(typeof(F))
    ω = unitless(uref(eom.f₀), 2π * F)
    RF_lin = Array{Complex{T}}(undef, eom.n_dof, eom.n_dof, n_fourier)
    RF_inp = Array{Complex{T}}(undef, eom.n_dof, size(eom.input, 2), n_fourier)
    RF_pnl = Array{Complex{T}}(undef, eom.n_dof, size(eom.p_nl, 2), n_fourier)
    for n in 1 : n_fourier
        RF = factorize(2.0 * (n * ω) * eom.K_lin[3] + im * eom.K_lin[2])
        RF_lin[:, :, n] .= im * (RF \ ((n * ω) ^ 2 * eom.K_lin[3] + im * n * ω * eom.K_lin[2] + eom.K_lin[1]))
        RF_inp[:, :, n] .= 2.0im * (RF \ eom.input)
        RF_pnl[:, :, n] .= -im * (RF \ eom.p_nl)
    end
    #DC = factorize(eom.K_lin[3])
    return RFSolver(
        eom, ω, n_fourier,
        RF_lin, RF_inp, RF_pnl,
        eom.M_lin,
        eom.M_p_nl
    )
end

RFSolver(:: Type{T}, circuit :: Circuit, F :: Frequency, f₀ :: Frequency, n_fourier :: Int = 1) where {T <: Real} = RFSolver(ClassicalEOM(T, circuit, f₀), F, n_fourier)
RFSolver(circuit :: Circuit, F :: Frequency, f₀ :: Frequency, n_fourier :: Int = 1) = RFSolver(ClassicalEOM(circuit, f₀), F, n_fourier)

Base.eltype(:: RFSolver{T}) where {T} = T
ndof(solver :: RFSolver) = solver.n_dof * (solver.n_fourier + 1) * 2

# dc current is -K₂ ẍ - K₁ ẋ + K₀ x - p_nl * sin(q_nl' * x - θ_nl)
# ẋ = v
# -K₂ v̇ - K₁ v + K₀ x - p_nl sin(q_nl' * x - θ_nl) = 0
# v̇ = -(K₂ \ K₁) * v + (K₂ \ K₀) x - (K₂ \ p_nl) sin(q_nl' * x - θ_nl)
# rf current i K' ẋ + K x - p_nl * sin(q_nl' * x - θ_nl) + 2 in
function rhs!(du, u, p :: Tuple{RFSolver, Function}, t)
    solver, v_i = p
    
    φ_dc = @view u[1 : solver.n_dof] 
    v_dc = @view u[solver.n_dof .+ (1 : solver.n_dof)]
    #φ_rf = reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof .+ (1 : 2 * solver.n_dof)])

    dφ_dc = @view du[1 : solver.n_dof] 
    dv_dc = @view du[solver.n_dof .+ (1 : solver.n_dof)]
    #dφ_rf = reinterpret(Complex{eltype(du)}, @view du[2 * solver.n_dof .+ (1 : 2 * solver.n_dof)])

    # linear dc part
    dφ_dc .= v_dc
    mul!(dv_dc, solver.DC_lin, @view u[1 : (2 * solver.n_dof)])

    # linear rf part
    for n in 1 : solver.n_fourier
        # n-th Fourier harmonic
        φ_rf = reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof * n .+ (1 : 2 * solver.n_dof)])
        dφ_rf = reinterpret(Complex{eltype(du)}, @view du[2 * solver.n_dof * n .+ (1 : 2 * solver.n_dof)])
        if n == 1
            # we have quasimonochromatic drive
            mul!(dφ_rf, (@view solver.RF_inp[:, :, n]), v_i(t))
            mul!(dφ_rf, (@view solver.RF_lin[:, :, n]), φ_rf, one(eltype(solver.RF_lin)), one(eltype(dφ_rf)))
        else
            mul!(dφ_rf, (@view solver.RF_lin[:, :, n]), φ_rf)
        end
    end

    # nonlinear part
    A = Vector{Complex{eltype(solver)}}(undef, 2 * solver.n_fourier + 1)
    for (i, (q, θ)) in enumerate(zip(eachcol(solver.q_nl), solver.θ_nl))
        # rf harmonics
        for n in 1 : solver.n_fourier
            A[1 + solver.n_fourier + n] = q' * reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof * n .+ (1 : 2 * solver.n_dof)])
            A[1 + solver.n_fourier - n] = conj(A[1 + solver.n_fourier + n])
        end
        # dc component
        A[1 + solver.n_fourier] = q' * φ_dc - θ
        m, B = inf_toeplitz_fun(-solver.n_fourier, A, sin)
         
        # dc part
        axpy!(real(B[-m + 1]), (@view solver.DC_p_nl[:, i]), dv_dc)
        # rf part
        for n in 1 : solver.n_fourier
            dφ_rf = reinterpret(Complex{eltype(du)}, @view du[2 * solver.n_dof * n .+ (1 : 2 * solver.n_dof)])
            axpy!(B[-m + 1 + n], (@view solver.RF_p_nl[:, i, n]), dφ_rf)
        end
    end
    return du
end

function jac!(jac, u, p :: Tuple{RFSolver, Function}, _)
    solver, _ = p
    fill!(jac, zero(eltype(jac)))
    n = solver.n_dof

    φ_dc = @view u[1 : n] 

    # linear dc part
    @simd for i in 1 : n
        jac[i, n + i] = one(eltype(jac))
    end
    jac[n .+ (1 : n), 1 : (2 * n)] .= solver.DC_lin

    # linear rf part
    for m in 1 : solver.n_fourier
        jac[2 * n * m .+ (1 : 2 * n), 2 * n * m .+ (1 : 2 : 2 * n)] .= reinterpret(eltype(solver), @view solver.RF_lin[:, :, m])
        jac[2 * n * m .+ (2 : 2 : 2 * n), 2 * n * m .+ (2 : 2 : 2 * n)] .= @view jac[2 * n * m .+ (1 : 2 : 2 * n), 2 * n * m .+ (1 : 2 : 2 * n)]
        jac[2 * n * m .+ (1 : 2 : 2 * n), 2 * n * m .+ (2 : 2 : 2 * n)] .-= @view jac[2 * n * m .+ (2 : 2 : 2 * n), 2 * n * m .+ (1 : 2 : 2 * n)]
    end

    # nonlinear part
    A = Vector{Complex{eltype(solver)}}(undef, 2 * solver.n_fourier + 1)
    for (i, (q, θ)) in enumerate(zip(eachcol(solver.q_nl), solver.θ_nl))
        # rf components
        for m in 1 : solver.n_fourier
            A[1 + solver.n_fourier + m] = q' * reinterpret(Complex{eltype(u)}, @view u[2 * solver.n_dof * m .+ (1 : 2 * solver.n_dof)])
            A[1 + solver.n_fourier - m] = conj(A[1 + solver.n_fourier + m])
        end
        # dc component
        A[1 + solver.n_fourier] = q' * φ_dc - θ
        l, B = inf_toeplitz_fun(-solver.n_fourier, A, cos)

        p_dc = @view solver.DC_p_nl[:, i]
        p_rf = @view solver.RF_p_nl[:, i, :]

        # dc part
        mul!(view(jac, n .+ (1 : n), 1 : n), p_dc, q', real(B[-l + 1]), one(eltype(jac)))
        for m in 1 : solver.n_fourier
            mul!(view(jac, n .+ (1 : n), (2 * n * m) .+ (1 : 2 : 2 * n)), p_dc, q', real(B[-l + 1 + m] + B[-l + 1 - m]), one(eltype(jac)))
            mul!(view(jac, n .+ (1 : n), (2 * n * m) .+ (2 : 2 : 2 * n)), p_dc, q', real(-im * B[-l + 1 + m] + im * B[-l + 1 - m]), one(eltype(jac)))
        end

        # rf part
        for m in 1 : solver.n_fourier
            # derivative wrt dc
            mul!(
                reinterpret(Complex{eltype(jac)}, view(jac, 2 * n * m .+ (1 : 2 * n), 1 : n)),
                (@view p_rf[:, m]), q', B[-l + 1 + m], one(Complex{eltype(jac)})
            )
            for m′ in 1 : solver.n_fourier
                B₊ = (1 <= (-l + 1 + m + m′) <= length(B)) ? B[-l + 1 + m + m′] : zero(eltype(B))
                B₋ = (1 <= (-l + 1 + m - m′) <= length(B)) ? B[-l + 1 + m - m′] : zero(eltype(B))
                mul!(
                    reinterpret(Complex{eltype(jac)}, view(jac, 2 * n * m .+ (1 : 2 * n), 2 * n * m′ .+ (1 : 2 : 2 * n))),
                    (@view p_rf[:, m]), q', B₊ + B₋, one(Complex{eltype(jac)})
                )
                mul!(
                    reinterpret(Complex{eltype(jac)}, view(jac, 2 * n * m .+ (1 : 2 * n), 2 * n * m′ .+ (2 : 2 : 2 * n))),
                    (@view p_rf[:, m]), q', -im * B₊ + im * B₋, one(Complex{eltype(jac)})
                )
            end
        end
    end
    return jac
end

"""
    steady_state(eom :: ClassicalEOM; n_newton :: Integer = 10)

Finds stationary solution of classical equations of motion `rf` using Newton method.
"""
steady_state(rf :: RFSolver{T}; n_newton :: Integer = 10) where {T} = [steady_state(rf.eom; n_newton); zeros(T, ndof(rf) - rf.n_dof)]

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
