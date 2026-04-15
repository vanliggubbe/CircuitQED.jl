abstract type AbstractSolver end
abstract type AbstractQuasiclassicalSolver <: AbstractSolver end
abstract type AbstractClassicalSolver <: AbstractQuasiclassicalSolver end

include("classical.jl")
include("rf.jl")
