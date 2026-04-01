const Mayhaps = Union{Nothing, T} where {T}

@inline uref(f :: Frequency) = (2π * f, 1.0u"ħ", 2.0u"q", 1.0u"k")
