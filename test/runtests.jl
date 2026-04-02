using CircuitQED
using Unitful
using LessUnits
using SciMLBase
using LinearAlgebra
using Test

@testset "CircuitQED.jl" begin
    unit_ref = (2π * 1.0u"GHz", 1.0u"ħ", 2.0u"q", 1.0u"k")
    @testset "Classical" begin
        circ = Circuit([
            Capacitor(:C_q, (ground(), :qubit), 200u"fF"),
            JoJunction(:J_q, (ground(), :qubit), 20u"GHz"),
            Capacitor(:C_c, (:qubit, :drive), 1u"fF"),
            Port(:drive, :drive, 50u"Ω")
        ])
            
        eom = ClassicalEOM(Float64, circ, 1u"GHz")
        @test eom.K_lin[1] ≈ zeros(2, 2)
        @test eom.K_lin[2] ≈ [0.0 0.0; 0.0 inv(unitless(unit_ref, 50u"Ω"))]
        @test eom.K_lin[3] ≈ [unitless(unit_ref, x) for x in [201u"fF" -1u"fF"; -1u"fF" 1u"fF"]]
        @test eom.p_nl[:, 1] * eom.q_nl[:, 1]' ≈ [20.0 0.0; 0.0 0.0]

        rf = RFSolver(eom, 10u"GHz")
        v_i(_) = [0.0]

        u = randn(ndof(rf))
        fun = ODEFunction(rf)
        jac_fd = zeros(eltype(rf), ndof(rf), ndof(rf))
        jac_ex = zeros(eltype(rf), ndof(rf), ndof(rf))
        @time fun.jac(jac_ex, u, (rf, v_i), 0.0)
        @time fun.jac(jac_ex, u, (rf, v_i), 0.0)
        ε = 1e-4
        for i in 1 : ndof(rf)
            v = copy(u)
            dv = similar(v)
   
            v[i] = u[i] + 0.5ε
            @time fun(dv, v, (rf, v_i), 0.0)
            jac_fd[:, i] .= dv
  
            v[i] = u[i] - 0.5ε
            @time fun(dv, v, (rf, v_i), 0.0)
            jac_fd[:, i] .-= dv
        end
        jac_fd ./= ε
        @test norm(jac_fd - jac_ex) < 1e-4
    end
end
