include("elements/elements.jl")

mutable struct Circuit
    nd_index :: Dict{Node, Int}
    el_index :: Dict{Tuple{Symbol, Symbol}, Int}
    nds :: Vector{Node}
    els :: Vector{Element}
end

Circuit() = Circuit(Dict(ground() => 1), Dict{Tuple{Symbol, Symbol}, Int}(), [ground()], Element[])

function Circuit(els)
    circuit = Circuit()
    add_elements!(circuit, els)
end

function _add_node!(circuit :: Circuit, node :: Node)
    if !haskey(circuit.nd_index, node)
        circuit.nd_index[node] = maximum(values(circuit.nd_index)) + 1
        push!(circuit.nds, node)
    end
    return circuit.nd_index[node]
end

function _compact_node_indices!(circuit :: Circuit)
    ids = sort!(collect(unique(values(circuit.nd_index))))
    renumber = Dict(id => i for (i, id) in enumerate(ids))
    for (node, id) in circuit.nd_index
        circuit.nd_index[node] = renumber[id]
    end

    nds = Vector{Node}(undef, length(ids))
    nds[1] = ground()
    assigned = falses(length(ids))
    assigned[1] = true
    for node in circuit.nds
        id = circuit.nd_index[node]
        if !assigned[id]
            nds[id] = node
            assigned[id] = true
        end
    end
    for (node, id) in circuit.nd_index
        if !assigned[id]
            nds[id] = node
            assigned[id] = true
        end
    end
    circuit.nds = nds
    return circuit
end

function add_element!(circuit :: Circuit, el :: Element)
    @argcheck !haskey(circuit.el_index, id(el)) "Element $(id(el)) already exists in the circuit"

    push!(circuit.els, el)
    circuit.el_index[id(el)] = length(circuit.els)

    for node in coordinates(el)
        _add_node!(circuit, node)
    end
    return circuit
end

function add_element!(circuit :: Circuit, el :: Short)
    @argcheck !haskey(circuit.el_index, id(el)) "Element $(id(el)) already exists in the circuit"

    push!(circuit.els, el)
    circuit.el_index[id(el)] = length(circuit.els)

    for node in coordinates(el)
        _add_node!(circuit, node)
    end

    ids = [circuit.nd_index[node] for node in coordinates(el)]
    target = any(==(circuit.nd_index[ground()]), ids) ? circuit.nd_index[ground()] : minimum(ids)
    losing = Set(filter(!=(target), ids))
    for (node, id) in circuit.nd_index
        if id in losing
            circuit.nd_index[node] = target
        end
    end
    _compact_node_indices!(circuit)
    return circuit
end

function add_elements!(circuit :: Circuit, els)
    for el in els
        add_element!(circuit, el)
    end
    return circuit
end

ndof(circuit :: Circuit) = maximum(values(circuit.nd_index)) - 1
@inline node_index(circuit, node :: Node) = circuit.nd_index[node] - 1
@inline node_index(circuit, node :: Node, default :: Int) = get(circuit.nd_index, node, default + 1) - 1

function _linear_idxs(circuit :: Circuit, el :: Element)
    idxs = Int[]
    idxs_el = Int[]
    for (i, node) in enumerate(coordinates(el))
        idx = node_index(circuit, node)
        if idx != 0
            push!(idxs, idx)
            push!(idxs_el, i)
        end
    end
    return idxs, idxs_el
end

# printing
function Base.show(io :: IO, circuit :: Circuit)
    println("Elements:")
    for el in filter(x -> !(x isa Port), circuit.els)
        println(io, "  ", el)
    end
    println("Ports:")
    for el in filter(x -> (x isa Port), circuit.els)
        println(io, "  ", el)
    end
end
