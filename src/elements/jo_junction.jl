struct JoJunction <: NonlinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    crit_current :: Ref{Mayhaps{Current}}
end

JoJunction(
    name,
    nodes :: Tuple{Any, Any},
    crit_current :: Mayhaps{<: Current} = nothing
) = JoJunction(Symbol(name), map(Node, nodes), Ref{Mayhaps{Current}}(crit_current))

JoJunction(
    name,
    nodes :: Tuple{Any, Any},
    j_energy :: Energy
) = JoJunction(Symbol(name), map(Node, nodes), Ref{Mayhaps{Current}}(uconvert(u"A", j_energy * 2π / 1.0u"Φ0")))

JoJunction(
    name,
    nodes :: Tuple{Any, Any},
    j_energy :: Frequency
) = JoJunction(Symbol(name), map(Node, nodes), Ref{Mayhaps{Current}}(uconvert(u"A", j_energy * 2π * 1.0u"h" / 1.0u"Φ0")))

dc_supercurrent(:: Type{<: JoJunction}) = (true, true)

nl_response(el :: JoJunction, f₀ :: Frequency) = (unitless(uref(f₀), el.crit_current[]) * [1.0, -1.0], [1.0, -1.0], [0.0])

Base.show(io :: IO, el :: JoJunction) = print(io, "JoJu $(el.name)($(el.nodes[1]), $(el.nodes[2])): Iᶜ = $(el.crit_current[])")
