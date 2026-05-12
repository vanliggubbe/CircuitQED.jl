struct CPWPiece <: LinearElement{3}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    impedance :: Ref{Mayhaps{ElectricalResistance}}
    hw_freq :: Ref{Mayhaps{Frequency}}
    n_aux :: Ref{Mayhaps{Integer}}
end

CPWPiece(
    name,
    nodes :: Tuple{Any, Any},
    impedance :: Mayhaps{<: ElectricalResistance} = nothing,
    hw_freq :: Mayhaps{<: Frequency} = nothing,
    n_aux :: Mayhaps{<: Integer} = nothing
) = CPWPiece(Symbol(name), map(Node, nodes), Ref{Mayhaps{ElectricalResistance}}(impedance), Ref{Mayhaps{Frequency}}(hw_freq), Ref{Mayhaps{Integer}}(n_aux))

dc_supercurrent(:: Type{<: CPWPiece}) = (false, true, true)

response(el :: CPWPiece, f₀ :: Frequency) = let Y = inv(unitless(uref(f₀), el.impedance[])), ω = unitless(uref(f₀), 2π * el.hw_freq[]), n = el.n_aux[]
    col1 = inv.(1 : n) * sqrt(2.0 * Y / π / ω)
    col2 = -col1 .* ((-1) .^ (1 : n))
    (
        [
            0.0         0.0             0.0             zeros(n)'
            0.0         (-ω * Y / π)    (ω * Y / π)     zeros(n)';
            0.0         (ω * Y / π)     (-ω * Y / π)    zeros(n)';
            zeros(n)    zeros(n)        zeros(n)        Diagonal(-(ω * (1 : n)) .^ 2)
        ],
        zeros(n + 3, n + 3),
        [
            0.0         0.0                 0.0                 zeros(n)'
            0.0         (π * Y / 3.0 / ω)   (π * Y / 6.0 / ω)   col1'
            0.0         (π * Y / 6.0 / ω)   (π * Y / 3.0 / ω)   col2'
            zeros(n)    col1                col2                1.0I(n)
        ]
    )
end

Base.show(io :: IO, el :: CPWPiece) = print(io, "CPW $(el.name)($(el.nodes[1]), $(el.nodes[2])): Z = $(el.impedance[]), f(λ/2) = $(el.hw_freq[]), # of modes = $(el.n_aux[])")
