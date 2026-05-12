struct Resistor <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    resistance :: Ref{Mayhaps{ElectricalResistance}}
    temperature :: Ref{Mayhaps{Temperature}}
end

Resistor(
    name,
    nodes :: Tuple{Any, Any},
    resistance :: Mayhaps{<: ElectricalResistance} = nothing,
    temperature :: Mayhaps{<: Temperature} = nothing
) = Resistor(Symbol(name), map(Node, nodes), Ref{Mayhaps{ElectricalResistance}}(resistance), Ref{Mayhaps{Temperature}}(temperature))

response(el :: Resistor, f₀ :: Frequency) = inv(unitless(uref(f₀), el.resistance[])) .* (zeros(2, 2), [1.0 -1.0; -1.0 1.0], zeros(2, 2))
