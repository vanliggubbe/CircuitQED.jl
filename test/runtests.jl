using CircuitQED
using Unitful
using LessUnits
using SciMLBase
using LinearAlgebra
using Test

function test_jacobian(eom; ε :: Real = 1e-4)
    v_i(_) = zeros(size(eom.dynamics.B, 2))
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
        @test ndof(eom) == size(eom.dynamics.A, 1)
        @test size(eom.dynamics.B, 2) == 1
        @test size(eom.output.A, 1) == 1
        @test eom.port_index[:drive] == 1

        u = test_jacobian(eom)
        @test length(output_voltage(eom, u, [0.0])) == 1
    end
    @testset "Classical reductions" begin
        massive = ClassicalEOM(Circuit([
            Capacitor(:C, (ground(), :n), 1u"fF"),
            Inductor(:L, (ground(), :n), 1u"nH"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massive) == 2
        @test massive.port_index[:p] == 1
        @test length(output_voltage(massive, test_jacobian(massive), [0.0])) == 1

        massless = ClassicalEOM(Circuit([
            Inductor(:L, (ground(), :n), 1u"nH"),
            Resistor(:R, (ground(), :n), 50u"Ω"),
            Port(:p, :n, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(massless) == 1
        @test massless.port_index[:p] == 1
        @test length(output_voltage(massless, test_jacobian(massless), [0.0])) == 1

        series = ClassicalEOM(Circuit([
            Inductor(:L1, (ground(), :mid), 1u"nH"),
            Inductor(:L2, (:mid, :out), 2u"nH"),
            Port(:p, :out, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(series) == 1
        @test series.port_index[:p] == 1
        @test length(output_voltage(series, test_jacobian(series), [0.0])) == 1

        capacitive_port = ClassicalEOM(Circuit([
            Capacitor(:C, (ground(), :node), 10u"fF"),
            Inductor(:L, (ground(), :node), 1u"nH"),
            Capacitor(:C_cpl, (:node, :input), 5u"fF"),
            Port(:drive, :input, 50u"Ω")
        ]), 1u"GHz")
        @test ndof(capacitive_port) == 3
        @test size(capacitive_port.coordinate_basis, 2) == 1
        @test capacitive_port.port_index[:drive] == 1
        @test length(output_voltage(capacitive_port, test_jacobian(capacitive_port), [0.0])) == 1
    end
end
