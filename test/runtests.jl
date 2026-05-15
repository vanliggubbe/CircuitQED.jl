using CircuitQED
using Unitful
using LessUnits
using SciMLBase
using LinearAlgebra
using Test

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

        v_i(_) = [0.0]
   
        u = randn(ndof(eom))
        fun = ODEFunction(eom)
        jac_fd = zeros(eltype(eom), ndof(eom), ndof(eom))
        jac_ex = zeros(eltype(eom), ndof(eom), ndof(eom))
        fun.jac(jac_ex, u, (eom, v_i), 0.0)
        ε = 1e-4
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
        #display(jac_ex)
        #display(jac_fd)
        @test norm(jac_fd - jac_ex) < 1e-4
        @test length(output_voltage(eom, u, [0.0])) == 1
    end
    @testset "RF" begin
        @test_throws ErrorException RFSolver(circ, 10u"GHz", 1u"GHz", 1)
    end
end
