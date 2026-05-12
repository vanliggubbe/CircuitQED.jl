struct Short{N} <: LinearElement{N}
    name :: Symbol
    nodes :: NTuple{N, Node}
end

Short(name, nodes :: Tuple) = Short(Symbol(name), map(Node, nodes))

dc_supercurrent(:: Type{<: Short{N}}) where {N} = tuple(fill(true, N)...)

Base.show(io :: IO, el :: Short) = print(io, "Short link $(el.name)($(join(el.nodes, ",")))")
