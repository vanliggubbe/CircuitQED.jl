using CircuitQED
using Unitful
using CairoMakie

CairoMakie.activate!()

rf2dc = Circuit([
    SNAIL(:S1, (ground(), :n1), 10u"μA", 3, 0.3, 0.35u"Φ0"),
    CPWPiece(:res, (:n1, :out), 50u"Ω", 10u"GHz", 3),
    Capacitor(:shunt, (:out, ground()), 2.5u"pF"),
    #CPWPiece(:filter, (:out, :open), 50u"Ω", 2 * 9.0u"GHz", 5),
    Port(:input, :n1, 50u"Ω"),
    Port(:output, :out, 50u"Ω")
])

eom = ClassicalEOM(Float64, rf2dc, 1u"GHz")

fs = LinRange(9, 11, 2000)
S = scattering_matrix(eom, fs * 1u"GHz")

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, fs, [abs2(s[1, 1]) for s in S])
lines!(ax1, fs, [abs2(s[2, 1]) for s in S])
ax2 = Axis(fig[2, 1])
lines!(ax2, fs, [angle(s[1, 1]) for s in S])
fig
