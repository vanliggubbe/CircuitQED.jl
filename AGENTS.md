# CircuitQED Agent Guidelines

## Development Commands
- Run all tests: `julia --project test/runtests.jl`
- Run tests with coverage: `julia --project -e 'using Pkg; Pkg.test(coverage=true)'`
- Format code: Not configured (JuliaFormatter recommended)
- Lint: Not configured (JuliaLint recommended)

## Project Structure
- Main package: `src/CircuitQED.jl`
- Solvers: `src/solvers/` (classical.jl, rf.jl)
- Elements: `src/elements/` (capacitor.jl, etc.)
- Utilities: `src/utils.jl`
- Tests: `test/runtests.jl`

## Testing
- Single test file: `test/runtests.jl` uses @testset blocks
- Test dependencies: Test, Unitful, LessUnits, SciMLBase, LinearAlgebra
- CI runs on Julia 1.10 and pre-release versions
- Tests create quantum circuits and verify Jacobians against finite differences

## Dependencies
- Core: LinearAlgebra, SparseArrays, SpecialFunctions
- Scientific: SciMLBase, FFTW, ArgCheck
- Units: Unitful, LessUnits
- Compatible with Julia 1.10+

## Notes
- Package uses Unitful.jl for physical quantities
- Follows Julia package conventions
- No pre-compiled artifacts or code generation