struct Capacitor <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    capacitance :: Ref{Mayhaps{Capacitance}}
end

Capacitor(name, nodes :: Tuple{Any, Any}, capacitance :: Mayhaps{<: Capacitance} = nothing) = Capacitor(Symbol(name), map(Node, nodes), Ref{Mayhaps{Capacitance}}(capacitance))

response(el :: Capacitor, f₀ :: Frequency) = unitless(uref(f₀), el.capacitance[]) .* (zeros(2, 2), zeros(2, 2), [1.0 -1.0; -1.0 1.0])

Base.show(io :: IO, el :: Capacitor) = print(io, "Capacitor $(el.name)($(el.nodes[1]), $(el.nodes[2])): C = $(el.capacitance[])")
