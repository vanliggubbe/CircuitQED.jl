using CircuitQED
using Unitful
using LessUnits
using SciMLBase
using LinearAlgebra
using Test
using NamedArrays
using OrdinaryDiffEq

function test_jacobian(eom; ε :: Real = 1e-4)
    v_i(_) = zeros(size(eom.dyn.B, 2))
    u = randn(ndof(eom))
    fun = ODEFunction(eom)
    jac_fd = zeros(eltype(eom), ndof(eom), ndof(eom))
    jac_ex = zeros(eltype(eom), ndof(eom), ndof(eom))
    fun.jac(jac_ex, u, (eom, v_i), 0.0)
    for i in 1 : ndof(eom)
        v = copy(u)
        dv = similar(v)

        v[i] = u[i] + 0.5ε
        fun(dv, v, (eom, v_i), 0.0)
        jac_fd[:, i] .= dv

        v[i] = u[i] - 0.5ε
        fun(dv, v, (eom, v_i), 0.0)
        jac_fd[:, i] .-= dv
    end
    jac_fd ./= ε
    @test norm(jac_fd - jac_ex) < ε
    return u
end

function relative_infnorm_error(numerical, exact)
    return maximum(abs, numerical .- exact) / maximum(abs, exact)
end

function test_linear_rcl_dynamics(C, L, R, v0, f, t_int, threshold)
    circ = Circuit([
        Capacitor(:C, (:node, ground()), C),
        Inductor(:L, (:node, ground()), L),
        Port(:R, :node, R)
    ])
    eom = ClassicalEOM(circ, 1u"GHz")
    units = eom.units
    Ω = unitless(units, 2π * f)
    γ = unitless(units, inv(R * C))
    ω0² = unitless(units, inv(L * C))
    v0_ul = unitless(units, v0)
    R_ul = unitless(units, R)

    A = [0.0 1.0; -ω0² -γ]
    b = [0.0, 2γ * v0_ul]
    z = (im * Ω * I - A) \ b
    x₀ = zeros(ndof(eom))

    v_i(t) = Dict(:R => v0 * sin(2π * f * t))
    prob = ODEProblem(eom, v_i, (0u"s", t_int), x₀)
    sol = solve(prob, Vern7(); reltol = 1e-10, abstol = 1e-12)

    ts = [t_int * (i - 1) / 299 for i in 1 : 300]
    v_exact = Float64[]
    i_exact = Float64[]
    for t in ts
        t_ul = unitless(units, t)
        x = imag.(z * exp(im * Ω * t_ul)) - exp(A * t_ul) * imag.(z)
        v_node = x[2]
        v_in = v0_ul * sin(Ω * t_ul)
        push!(v_exact, v_node - v_in)
        push!(i_exact, (v_node + v_in) / R_ul)
    end

    v_ode = [unitless(units, output_voltage(sol, :R, t)) for t in ts]
    i_ode = [unitless(units, output_current(sol, :R, t)) for t in ts]

    @test relative_infnorm_error(v_ode, v_exact) < threshold
    @test relative_infnorm_error(i_ode, i_exact) < threshold
end

@testset "CircuitQED.jl" begin
    unit_ref = (2π * 1.0u"GHz", 1.0u"ħ", 2.0u"q", 1.0u"k")
    circ = Circuit([
        Capacitor(:C_q, (ground(), :qubit), 200u"fF"),
        JoJunction(:J_q, (ground(), :qubit), 20u"GHz"),
        Capacitor(:C_c, (:qubit, :drive), 1u"fF"),
        Port(:drive, :drive, 50u"Ω")
    ])

    @testset "Classical" begin
        eom = ClassicalEOM(circ, 1u"GHz")
        @test ndof(eom) == size(eom.dyn.A, 1)
        @test size(eom.dyn.B, 2) == 1
        @test eom.circuit === circ
        @test eom.circuit.port_index[:drive] == 1
        @test eom.circuit.port_list[1] isa Port

        u = test_jacobian(eom)
    end
    @testset "Classical reductions" begin
        massive = ClassicalEOM(Circuit([
            Capacitor(:C, (ground(), :n), 1u"fF"),
            Inductor(:L, (ground(), :n), 1u"nH"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massive) == 2
        @test massive.circuit.port_index[:p] == 1
        @test massive.circuit.port_list[1].name == :p

        massless = ClassicalEOM(Circuit([
            Inductor(:L, (ground(), :n), 1u"nH"),
            Resistor(:R, (ground(), :n), 50u"Ω"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massless) == 1
        @test massless.circuit.port_index[:p] == 1

        series = Circuit([
            Inductor(:L1, (ground(), :mid), 1u"nH"),
            Inductor(:L2, (:mid, :out), 2u"nH"),
            Port(:p, :out, 50u"Ω")
        ])
        @test_throws ErrorException ClassicalEOM(series, 1u"GHz")

        capacitive_port = ClassicalEOM(Circuit([
            Capacitor(:C, (ground(), :node), 10u"fF"),
            Inductor(:L, (ground(), :node), 1u"nH"),
            Capacitor(:C_cpl, (:node, :input), 5u"fF"),
            Port(:drive, :input, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(capacitive_port) == size(capacitive_port.dyn.A, 1)
        @test capacitive_port.circuit.port_index[:drive] == 1

        u = test_jacobian(capacitive_port)
        v_in = [0.25]
        v_out = CircuitQED._output_voltage(capacitive_port, u, v_in)
        i_out = CircuitQED._output_current(capacitive_port, u, v_in)
        Y = inv(unitless(capacitive_port.units, capacitive_port.circuit.port_list[1].impedance[]))
        @test i_out ≈ Y .* (v_out .+ 2 .* v_in)
    end
    @testset "S-matrix" begin
        # a linear circuit with both damped mode and cyclic phase
        C = 50u"fF"
        L = 1u"nH"
        L_cpl = 0.05u"nH"
        C_cpl = 1u"fF"
        Z = 50u"Ω"
        Ic = 20u"nA"

        circ = Circuit([
            Capacitor(:C, (:node, ground()), C),
            Inductor(:L, (:node, :cpl_1), L),
            Inductor(:L_cpl, (:cpl_1, ground()), L_cpl),
            Capacitor(:C_cpl, (:node, :cpl_2), C_cpl),
            Port(:port_1, :cpl_1, Z),
            Port(:port_2, :cpl_2, Z)
        ])
        eom = ClassicalEOM(circ, 1u"GHz")

        function S_matrix_1(f :: Unitful.Frequency)
            ω = unitless(eom.units, 2π * f)
            right = -[0.0 0.0; 2.0 0.0; 0.0 2.0] / unitless(eom.units, Z)
            left = [0.0 (-im * ω) 0.0; 0.0 0.0 (-im * ω)]
            M0 = [
                unitless(eom.units, -inv(L))        unitless(eom.units, inv(L))                 0.0;
                unitless(eom.units, inv(L))         unitless(eom.units, -inv(L) - inv(L_cpl))   0.0;
                0.0                                 0.0                                         0.0
            ]
            M1 = [
                0.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
            ] / unitless(eom.units, Z)
            M2 = [
                unitless(eom.units, C + C_cpl)  0.0 unitless(eom.units, -C_cpl);
                0.0                             0.0 0.0;
                unitless(eom.units, -C_cpl)     0.0 unitless(eom.units, C_cpl)
            ]
            return NamedArray(left * ((M0 + im * ω * M1 + ω ^ 2 * M2) \ right) - I, ([:port_1, :port_2], [:port_1, :port_2]))
        end
        fs = LinRange(-10u"GHz", 10u"GHz", 20)
        Ss = scattering_matrix(eom, fs)
        # test S matrix
        @test all([S ≈ S_matrix_1(f) for (S, f) in zip(Ss, fs)])

        circ2 = Circuit([
            Capacitor(:C, (:node, ground()), C),
            Inductor(:L, (:node, :cpl_1), L),
            JoJunction(:J, (:node, ground()), Ic),
            Inductor(:L_cpl, (:cpl_1, ground()), L_cpl),
            Capacitor(:C_cpl, (:node, :cpl_2), C_cpl),
            Port(:port_1, :cpl_1, Z),
            Port(:port_2, :cpl_2, Z)
        ])
        eom2 = ClassicalEOM(circ2, 1u"GHz")

        function S_matrix_2(f :: Unitful.Frequency)
            ω = unitless(eom.units, 2π * f)
            right = -[0.0 0.0; 2.0 0.0; 0.0 2.0] / unitless(eom.units, Z)
            left = [0.0 (-im * ω) 0.0; 0.0 0.0 (-im * ω)]
            M0 = [
                unitless(eom.units, -inv(L))        unitless(eom.units, inv(L))                 0.0;
                unitless(eom.units, inv(L))         unitless(eom.units, -inv(L) - inv(L_cpl))   0.0;
                0.0                                 0.0                                         0.0
            ]
            M0[1, 1] -= unitless(eom.units, Ic)
            M1 = [
                0.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
            ] / unitless(eom.units, Z)
            M2 = [
                unitless(eom.units, C + C_cpl)  0.0 unitless(eom.units, -C_cpl);
                0.0                             0.0 0.0;
                unitless(eom.units, -C_cpl)     0.0 unitless(eom.units, C_cpl)
            ]
            return NamedArray(left * ((M0 + im * ω * M1 + ω ^ 2 * M2) \ right) - I, ([:port_1, :port_2], [:port_1, :port_2]))
        end
        fs = LinRange(-10u"GHz", 10u"GHz", 20)
        Ss = scattering_matrix(eom2, fs)
        # test S matrix
        @test all([S ≈ S_matrix_2(f) for (S, f) in zip(Ss, fs)])
    end
    @testset "Linear time dynamics" begin
        test_linear_rcl_dynamics(50u"fF", 1u"nH", 50u"Ω", 1u"μV", 2u"GHz", 20u"ns", 1e-5)
        test_linear_rcl_dynamics(50u"fF", 1u"nH", 50u"kΩ", 1u"μV", 20u"GHz", 20u"ns", 1e-5)
        test_linear_rcl_dynamics(80u"fF", 2u"nH", 1000u"Ω", 2u"μV", 5u"GHz", 30u"ns", 1e-5)
        test_linear_rcl_dynamics(200u"fF", 0.2u"nH", 30u"kΩ", 2u"μV", 5u"GHz", 30u"ns", 1e-5)
    end
end
