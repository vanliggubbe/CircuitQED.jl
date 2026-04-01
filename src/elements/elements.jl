const Node = Symbol
abstract type Element{N} end
abstract type LinearElement{N} <: Element{N} end
abstract type NonlinearElement{N} <: Element{N} end

ground() = :ground

struct Capacitor <: LinearElement{2} 
    name :: Symbol
    nodes :: Tuple{Node, Node}
    capacitance :: Ref{Mayhaps{Capacitance}}
end

struct Inductor <: LinearElement{2} 
    name :: Symbol
    nodes :: Tuple{Node, Node}
    inductance :: Ref{Mayhaps{Inductance}}
end

struct Resistor <: LinearElement{2} 
    name :: Symbol
    nodes :: Tuple{Node, Node}
    resistance :: Ref{Mayhaps{ElectricalResistance}}
    temperature :: Ref{Mayhaps{Temperature}}
end

struct Port <: LinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node}
    impedance :: Ref{Mayhaps{ElectricalResistance}}
    temperature :: Ref{Mayhaps{Temperature}}
end

struct CPWPiece <: LinearElement{3}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    impedance :: Ref{Mayhaps{ElectricalResistance}}
    hw_freq :: Ref{Mayhaps{Frequency}}
    n_aux :: Ref{Mayhaps{Integer}}
end

struct JoJunction <: NonlinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    crit_current :: Ref{Mayhaps{Current}}
end

struct SNAIL <: NonlinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    crit_current :: Ref{Mayhaps{Current}}
    n :: Ref{Mayhaps{Int}}
    asymmetry :: Ref{Mayhaps{Real}}
    flux :: Ref{Mayhaps{MagneticFlux}}
end

# constructors 
Capacitor(name :: Symbol, nodes :: Tuple{Node, Node}, capacitance :: Mayhaps{<: Capacitance} = nothing) = Capacitor(name, nodes, Ref{Mayhaps{Capacitance}}(capacitance))

Inductor(name :: Symbol, nodes :: Tuple{Node, Node}, inductance :: Mayhaps{<: Inductance} = nothing) = Inductor(name, nodes, Ref{Mayhaps{Inductance}}(inductance))

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

JoJunction(
    name :: Symbol,
    nodes :: Tuple{Node, Node},
    crit_current :: Mayhaps{<: Current} = nothing
) = JoJunction(name, nodes, Ref{Mayhaps{Current}}(crit_current))

JoJunction(
    name :: Symbol,
    nodes :: Tuple{Node, Node},
    j_energy :: Energy
) = JoJunction(name, nodes, Ref{Mayhaps{Current}}(uconvert(u"A", j_energy * 2π / 1.0u"Φ0")))

JoJunction(
    name :: Symbol,
    nodes :: Tuple{Node, Node},
    j_energy :: Frequency
) = JoJunction(name, nodes, Ref{Mayhaps{Current}}(uconvert(u"A", j_energy * 2π * 1.0u"h" / 1.0u"Φ0")))

SNAIL(
    name :: Symbol,
    nodes :: Tuple{Node, Node},
    crit_current :: Mayhaps{<: Current} = nothing,
    n :: Mayhaps{Int} = nothing,
    asymmetry :: Mayhaps{<: Real} = nothing,
    flux :: Mayhaps{<: MagneticFlux} = nothing
) = SNAIL(name, nodes, Ref{Mayhaps{Current}}(crit_current), Ref{Mayhaps{Int}}(n), Ref{Mayhaps{Real}}(asymmetry), Ref{Mayhaps{MagneticFlux}}(flux))


is_lossy(:: Type{T}) where {T <: Element} = hasfield(T, :temperature)
is_lossless(:: Type{T}) where {T <: Element} = !is_lossless(T)
is_linear(:: Type{<: LinearElement}) = true
is_linear(:: Type{<: NonlinearElement}) = false
is_nonlinear(:: Type{<: LinearElement}) = false
is_nonlinear(:: Type{<: NonlinearElement}) = true

dc_supercurrent(:: Type{<: Element{2}}) = (false, false)
dc_supercurrent(:: Type{<: Inductor}) = (true, true)
dc_supercurrent(:: Type{<: JoJunction}) = (true, true)
dc_supercurrent(:: Type{<: SNAIL}) = (true, true)
dc_supercurrent(:: Type{<: CPWPiece}) = (false, true, true)

for fun in [:is_lossy, :is_lossless, :is_linear, :is_nonlinear, :dc_supercurrent]
    @eval $(fun)(:: T) where {T <: Element} = $(fun)(T)
end

id(el :: T) where {T <: Element} = (Symbol(T), el.name)

nodes(el :: Element{N}) where {N} = (length(el.nodes) == N ? el.nodes : (ground(), el.nodes...))

# response functions
response(el :: Capacitor,   f₀ :: Frequency) = unitless(uref(f₀), el.capacitance[])     .* (zeros(2, 2), zeros(2, 2), [1.0 -1.0; -1.0 1.0])
response(el :: Inductor,    f₀ :: Frequency) = inv(unitless(uref(f₀), el.inductance[])) .* ([-1.0 1.0; 1.0 -1.0], zeros(2, 2), zeros(2, 2))
response(el :: Resistor,    f₀ :: Frequency) = inv(unitless(uref(f₀), el.resistance[])) .* (zeros(2, 2), [1.0 -1.0; -1.0 1.0], zeros(2, 2))
response(el :: Port,        f₀ :: Frequency) = inv(unitless(uref(f₀), el.impedance[]))  .* (zeros(2, 2), [0.0 0.0; 0.0 1.0], zeros(2, 2))
response(el :: CPWPiece,    f₀ :: Frequency) = let Y = inv(unitless(uref(f₀), el.impedance[])), ω = unitless(uref(f₀), 2π * el.hw_freq[]), n = el.n_aux[]
    col1 = inv.(1 : n) * sqrt(2.0 * Y / π / ω)
    col2 = col1 .* ((-1) .^ (1 : n))
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

# nonlinear elements response
nl_response(el :: JoJunction,   f₀ :: Frequency) = (unitless(uref(f₀), el.crit_current[]) * [1.0, -1.0], [1.0, -1.0], [0.0])
nl_response(el :: SNAIL,        f₀ :: Frequency) = (
    unitless(uref(f₀), el.crit_current[]) * [
        1.0     el.asymmetry[];
        -1.0    -el.asymmetry[]
    ], [
        1.0     inv(el.n[]);
        -1.0    -inv(el.n[])
    ], [0.0, unitless(uref(f₀), el.flux[]) / el.n[]]
)

#printing
Base.show(io :: IO, el :: Capacitor)    = print(io, "Capacitor $(el.name)($(el.nodes[1]), $(el.nodes[2])): C = $(el.capacitance[])")
Base.show(io :: IO, el :: Inductor)     = print(io, "Inductor $(el.name)($(el.nodes[1]), $(el.nodes[2])): L = $(el.inductance[])")
Base.show(io :: IO, el :: Port)         = print(io, "Port $(el.name)($(el.nodes[1])): Z = $(el.impedance[]), T = $(el.temperature[])")
Base.show(io :: IO, el :: CPWPiece)     = print(io, "CPW $(el.name)($(el.nodes[1]), $(el.nodes[2])): Z = $(el.impedance[]), f(λ/2) = $(el.hw_freq[]), # of modes = $(el.n_aux[])")
Base.show(io :: IO, el :: JoJunction)   = print(io, "JoJu $(el.name)($(el.nodes[1]), $(el.nodes[2])): Iᶜ = $(el.crit_current[])")
Base.show(io :: IO, el :: SNAIL)        = print(io, "SNAIL $(el.name)($(el.nodes[1]), $(el.nodes[2])): Iᶜ = $(el.crit_current[]), α = $(el.asymmetry[]), n = $(el.n[])")
