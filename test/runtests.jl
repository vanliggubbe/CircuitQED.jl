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
        @test eom.K_lin[1] ≈ zeros(2, 2)
        @test eom.K_lin[2] ≈ [0.0 0.0; 0.0 inv(unitless(unit_ref, 50u"Ω"))]
        @test eom.K_lin[3] ≈ [unitless(unit_ref, x) for x in [201u"fF" -1u"fF"; -1u"fF" 1u"fF"]]
        @test eom.p_nl[:, 1] * eom.q_nl[:, 1]' ≈ [20.0 0.0; 0.0 0.0]

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
    end
    @testset "RF" begin
        for n_fourier in 1 : 2
            rf = RFSolver(circ, 10u"GHz", 1u"GHz", n_fourier)
            v_i(_) = [0.0]
    
            u = randn(ndof(rf))
            fun = ODEFunction(rf)
            jac_fd = zeros(eltype(rf), ndof(rf), ndof(rf))
            jac_ex = zeros(eltype(rf), ndof(rf), ndof(rf))
            fun.jac(jac_ex, u, (rf, v_i), 0.0)
            ε = 1e-4
            for i in 1 : ndof(rf)
                v = copy(u)
                dv = similar(v)
       
                v[i] = u[i] + 0.5ε
                fun(dv, v, (rf, v_i), 0.0)
                jac_fd[:, i] .= dv
      
                v[i] = u[i] - 0.5ε
                fun(dv, v, (rf, v_i), 0.0)
                jac_fd[:, i] .-= dv
            end
            jac_fd ./= ε
            @test norm(jac_fd - jac_ex) < 1e-4
        end
    end
end
