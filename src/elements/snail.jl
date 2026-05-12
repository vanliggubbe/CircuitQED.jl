struct SNAIL <: NonlinearElement{2}
    name :: Symbol
    nodes :: Tuple{Node, Node}
    crit_current :: Ref{Mayhaps{Current}}
    n :: Ref{Mayhaps{Int}}
    asymmetry :: Ref{Mayhaps{Real}}
    flux :: Ref{Mayhaps{MagneticFlux}}
end

SNAIL(
    name :: Symbol,
    nodes :: Tuple{Node, Node},
    crit_current :: Mayhaps{<: Current} = nothing,
    n :: Mayhaps{Int} = nothing,
    asymmetry :: Mayhaps{<: Real} = nothing,
    flux :: Mayhaps{<: MagneticFlux} = nothing
) = SNAIL(name, nodes, Ref{Mayhaps{Current}}(crit_current), Ref{Mayhaps{Int}}(n), Ref{Mayhaps{Real}}(asymmetry), Ref{Mayhaps{MagneticFlux}}(flux))

dc_supercurrent(:: Type{<: SNAIL}) = (true, true)

nl_response(el :: SNAIL, f₀ :: Frequency) = (
    unitless(uref(f₀), el.crit_current[]) * [
        1     el.asymmetry[];
        -1    -el.asymmetry[]
    ], [
        inv(el.n[])     1;
        -inv(el.n[])    -1;
    ], [unitless(uref(f₀), el.flux[]) / el.n[], 0]
)

Base.show(io :: IO, el :: SNAIL) = print(io, "SNAIL $(el.name)($(el.nodes[1]), $(el.nodes[2])): Iᶜ = $(el.crit_current[]), α = $(el.asymmetry[]), n = $(el.n[])")
