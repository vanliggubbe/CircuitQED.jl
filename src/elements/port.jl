struct Port <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node}
    impedance :: Ref{Mayhaps{ElectricalResistance}}
    temperature :: Ref{Mayhaps{Temperature}}
end

Port(
    name :: Symbol,
    nodes :: Tuple{Node},
    impedance :: Mayhaps{<: ElectricalResistance} = nothing,
    temperature :: Mayhaps{<: Temperature} = nothing
) = Port(name, nodes, Ref{Mayhaps{ElectricalResistance}}(impedance), Ref{Mayhaps{Temperature}}(temperature))

Port(
    name :: Symbol,
    node :: Node,
    impedance :: Mayhaps{<: ElectricalResistance} = nothing,
    temperature :: Mayhaps{<: Temperature} = nothing
) = Port(name, (node,), Ref{Mayhaps{ElectricalResistance}}(impedance), Ref{Mayhaps{Temperature}}(temperature))

response(el :: Port, f₀ :: Frequency) = inv(unitless(uref(f₀), el.impedance[])) .* (zeros(2, 2), [0.0 0.0; 0.0 1.0], zeros(2, 2))

Base.show(io :: IO, el :: Port) = print(io, "Port $(el.name)($(el.nodes[1])): Z = $(el.impedance[]), T = $(el.temperature[])")
