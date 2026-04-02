using CircuitQED
using Unitful
using LessUnits
using CairoMakie
using OrdinaryDiffEq

CairoMakie.activate!()

rf2dc = Circuit([
    SNAIL(:S1, (ground(), :n3), 10u"μA", 3, 0.3, 0.4u"Φ0"),
    Capacitor(:C3, (ground(), :n3), 1u"fF"),
    SNAIL(:S3, (:n3, :n2), 10u"μA", 3, 0.3, 0.4u"Φ0"),
    Capacitor(:C2, (ground(), :n2), 1u"fF"),
    SNAIL(:S2, (:n2, :n1), 10u"μA", 3, 0.3, 0.4u"Φ0"),
    CPWPiece(:res, (:n1, :out), 50u"Ω", 10u"GHz", 3),
    Capacitor(:shunt, (:out, ground()), 4.0u"pF"),
    #CPWPiece(:filter, (:out, :open), 50u"Ω", 2 * 9.0u"GHz", 5),
    Port(:input, :n1, 50u"Ω"),
    Port(:output, :out, 50u"Ω")
])

eom = ClassicalEOM(Float64, rf2dc, 1u"GHz")

fs = LinRange(9, 11, 2000)
S = scattering_matrix(eom, fs * 1u"GHz")

rf = RFSolver(eom, fs[argmax([abs(s[2, 1]) for s in S])] * 1.0u"GHz")
prob = ODEProblem(rf,  t -> ([exp(-(t / 20u"ns") ^ 2), 0.0] * 120u"μV"), (-200u"ns", 200u"ns"))
#fun = ODEFunction(rf)
#display(fun(zeros(4 * eom.n_dof), steady_state(rf), (rf, t -> [0.0, 0.0]), 0.0))
sol = solve(prob, Rodas5P(); reltol = 1e-9, abstol = 1e-9)

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, fs, [abs2(s[1, 1]) for s in S])
lines!(ax1, fs, [abs2(s[2, 1]) for s in S])
ax2 = Axis(fig[2, 1])
lines!(ax2, sol.t / unitless(CircuitQED.uref(1u"GHz"), 1u"ns"), [u[eom.n_dof + 2] for u in sol.u] / unitless(CircuitQED.uref(1u"GHz"), 1u"μV"))
fig
