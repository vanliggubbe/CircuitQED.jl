const Node = Symbol
abstract type Element{N} end
abstract type LinearElement{N} <: Element{N} end
abstract type NonlinearElement{N} <: Element{N} end

ground() = :ground

is_lossy(:: Type{T}) where {T <: Element} = hasfield(T, :temperature)
is_lossless(:: Type{T}) where {T <: Element} = !is_lossy(T)
is_linear(:: Type{<: LinearElement}) = true
is_linear(:: Type{<: NonlinearElement}) = false
is_nonlinear(:: Type{<: LinearElement}) = false
is_nonlinear(:: Type{<: NonlinearElement}) = true

dc_supercurrent(:: Type{<: Element{2}}) = (false, false)

for fun in [:is_lossy, :is_lossless, :is_linear, :is_nonlinear, :dc_supercurrent]
    @eval $(fun)(:: T) where {T <: Element} = $(fun)(T)
end

id(el :: T) where {T <: Element} = (Symbol(T), el.name)

terminals(el :: Element) = el.nodes
coordinates(el :: Element{N}) where {N} = (length(el.nodes) == N ? el.nodes : (ground(), el.nodes...))
nodes(el :: Element) = coordinates(el)

has_inertia(el :: LinearElement, f₀ :: Frequency) = !all(iszero, response(el, f₀)[3])
has_damping(el :: LinearElement, f₀ :: Frequency) = !all(iszero, response(el, f₀)[2])
has_potential(el :: LinearElement, f₀ :: Frequency) = !all(iszero, response(el, f₀)[1])
has_inertia(el :: NonlinearElement, f₀ :: Frequency) = false
has_damping(el :: NonlinearElement, f₀ :: Frequency) = false
has_potential(el :: NonlinearElement, f₀ :: Frequency) = true

include("capacitor.jl")
include("inductor.jl")
include("resistor.jl")
include("port.jl")
include("cpw_piece.jl")
include("jo_junction.jl")
include("snail.jl")
include("short.jl")
