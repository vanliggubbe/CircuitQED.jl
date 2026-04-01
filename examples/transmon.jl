using CircuitQED
using Unitful

circ = Circuit([
    Capacitor(:C_q, (ground(), :qubit), 200u"fF"),
    JoJunction(:J_q, (ground(), :qubit), 20u"GHz"),
    Capacitor(:C_c, (:qubit, :drive), 1u"fF"),
    Port(:drive, :drive, 50u"Ω")
])
        
eom = ClassicalEOM(Float64, circ, 1u"GHz")
rf = RFSolver(eom, 12u"GHz")

v_i(_) = [0.0]

