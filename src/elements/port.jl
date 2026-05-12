struct Port <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node}
    impedance :: Ref{Mayhaps{ElectricalResistance}}
    temperature :: Ref{Mayhaps{Temperature}}
end

Port(
    name,
    nodes :: Tuple{Any},
    impedance :: Mayhaps{<: ElectricalResistance} = nothing,
    temperature :: Mayhaps{<: Temperature} = nothing
) = Port(Symbol(name), map(Node, nodes), Ref{Mayhaps{ElectricalResistance}}(impedance), Ref{Mayhaps{Temperature}}(temperature))

Port(
    name,
    node,
    impedance :: Mayhaps{<: ElectricalResistance} = nothing,
    temperature :: Mayhaps{<: Temperature} = nothing
) = Port(Symbol(name), (Node(node),), Ref{Mayhaps{ElectricalResistance}}(impedance), Ref{Mayhaps{Temperature}}(temperature))

response(el :: Port, f₀ :: Frequency) = inv(unitless(uref(f₀), el.impedance[])) .* (zeros(2, 2), [0.0 0.0; 0.0 1.0], zeros(2, 2))

Base.show(io :: IO, el :: Port) = print(io, "Port $(el.name)($(el.nodes[1])): Z = $(el.impedance[]), T = $(el.temperature[])")
