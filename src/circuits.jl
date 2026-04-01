include("elements/elements.jl")

mutable struct Circuit
    adj :: SparseMatrixCSC{Bool, Int}
    nd_index :: Dict{Node, Int}
    el_index :: Dict{Tuple{Symbol, Symbol}, Int}
    nds :: Vector{Node}
    els :: Vector{Element}
end

Circuit() = Circuit(sparse(Matrix{Bool}(undef, 0, 1)), Dict(ground() => 1), Dict{Tuple{Symbol, Symbol}, Int}(), [ground()], Element[])

function Circuit(els)
    circuit = Circuit()
    add_elements!(circuit, els)
end

function add_element!(circuit :: Circuit, el :: Element)
    if haskey(circuit.el_index, id(el))
        error("Element is there")
    end

    push!(circuit.els, el)
    circuit.el_index[id(el)] = length(circuit.els)

    for node in filter(x -> !haskey(circuit.nd_index, x), nodes(el))
        push!(circuit.nds, node)
        circuit.nd_index[node] = length(circuit.nds)
    end
    
    circuit.adj = sparse(findnz(circuit.adj)..., length(circuit.els), length(circuit.nds))
    for node in nodes(el)
        circuit.adj[length(circuit.els), circuit.nd_index[node]] = true
    end
    return circuit
end

function add_elements!(circuit :: Circuit, els)
    for el in els
        add_element!(circuit, el)
    end
    return circuit
end

ndof(circuit :: Circuit) = (length(circuit.nd_index) - 1 + sum((hasfield(typeof(el), :n_aux) ? el.n_aux[] : 0) for el in circuit.els))
@inline node_index(circuit, node :: Node) = circuit.nd_index[node] - 1
@inline node_index(circuit, node :: Node, default :: Int) = get(circuit.nd_index, node, default + 1) - 1

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
