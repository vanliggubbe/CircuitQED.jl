using CircuitQED
using Unitful
using LessUnits
using SciMLBase
using LinearAlgebra
using Test
using NamedArrays

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
        @test eom.port_index[:drive] == 1

        u = test_jacobian(eom)
    end
    @testset "Classical reductions" begin
        massive = ClassicalEOM(Circuit([
            Capacitor(:C, (ground(), :n), 1u"fF"),
            Inductor(:L, (ground(), :n), 1u"nH"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massive) == 2
        @test massive.port_index[:p] == 1

        massless = ClassicalEOM(Circuit([
            Inductor(:L, (ground(), :n), 1u"nH"),
            Resistor(:R, (ground(), :n), 50u"Ω"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massless) == 1
        @test massless.port_index[:p] == 1

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
        @test capacitive_port.port_index[:drive] == 1
    end
    @testset "Linear circuits" begin
        # a linear circuit with both damped mode and cyclic phase
        C = 50u"fF"
        L = 1u"nH"
        L_cpl = 0.05u"nH"
        C_cpl = 1u"fF"
        Z = 50u"Ω"
        circ = Circuit([
            Capacitor(:C, (:node, ground()), C),
            Inductor(:L, (:node, :cpl_1), L),
            Inductor(:L_cpl, (:cpl_1, ground()), L_cpl),
            Capacitor(:C_cpl, (:node, :cpl_2), C_cpl),
            Port(:port_1, :cpl_1, Z),
            Port(:port_2, :cpl_2, Z)
        ])
        eom = ClassicalEOM(circ, 1u"GHz")

        function S_matrix(f :: Unitful.Frequency)
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
        @test all([S ≈ S_matrix(f) for (S, f) in zip(Ss, fs)])
    end
end
