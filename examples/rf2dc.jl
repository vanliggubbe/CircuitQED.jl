using CircuitQED
using Unitful

rf2dc = Circuit([
    SNAIL(:S1, (ground(), :n1), 10u"μA", 3, 0.3, 0.35u"Φ0"),
    CPWPiece(:res, (:n1, :out), 50u"Ω", 10u"GHz", 3),
    Capacitor(:shunt, (:out, ground()), 1.5u"pF"),
    Port(:input, :n1, 50u"Ω"),
    Port(:output, :out, 50u"Ω")
])

eom = ClassicalEOM(Float64, rf2dc, 1u"GHz")
