struct Inductor <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    inductance :: Ref{Mayhaps{Inductance}}
end

Inductor(name, nodes :: Tuple{Any, Any}, inductance :: Mayhaps{<: Inductance} = nothing) = Inductor(Symbol(name), map(Node, nodes), Ref{Mayhaps{Inductance}}(inductance))

dc_supercurrent(:: Type{<: Inductor}) = (true, true)

response(el :: Inductor, f₀ :: Frequency) = inv(unitless(uref(f₀), el.inductance[])) .* ([-1.0 1.0; 1.0 -1.0], zeros(2, 2), zeros(2, 2))

Base.show(io :: IO, el :: Inductor) = print(io, "Inductor $(el.name)($(el.nodes[1]), $(el.nodes[2])): L = $(el.inductance[])")
