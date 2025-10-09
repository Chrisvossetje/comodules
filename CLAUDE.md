# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust library (`comodules`) for calculating Ext for coalgebras and comodules using high-level homological algebra. The library works over finite fields and (graded) univariate polynomial rings over a field. The codebase includes both a library and a binary executable for computational algebra operations.

## Build and Development Commands

### Building and Running
- `cargo build` - Build the library
- `cargo run` - Run the main binary (processes coalgebra examples)
- `cargo build --release` - Build optimized release version
- `cargo build --profile=profile` - Build with debug symbols for profiling

### Testing
- `cargo test` - Run all tests
- `cargo test <module_name>` - Run tests for specific module (e.g., `cargo test kcoalgebra_tests`)
- `cargo test -- --nocapture` - Run tests with output visible

### Library vs Binary
- The library is named `comodules` (src/lib.rs)
- The binary is named `comodule` (src/bin/main.rs)

## Core Architecture

### Module Structure
The codebase is organized into several main modules:

1. **`comodule/`** - Core algebraic structures
   - `kcoalgebra` - Coalgebra implementations over finite fields
   - `kcomodule` - Comodule implementations with k-linear structure
   - `rcomodule` - R-linear comodule implementations
   - `kmorphism`/`rmorphism` - Morphisms between comodules
   - `tensor` - Tensor product operations
   - `parsers` - Text file parsing for algebraic structures
   - `traits` - Core trait definitions (`Comodule`, `ComoduleMorphism`)

2. **`linalg/`** - Linear algebra foundations
   - `field` - Finite field implementations (primarily F2)
   - `flat_matrix`/`row_matrix` - Matrix implementations
   - `graded`/`grading` - Graded structure support
   - `ring`/`module` - Ring and module theory (supports polynomial rings over fields)
   - `groups` - Group theory components

3. **`polynomial/`** - Polynomial algebra
   - `buchberger` - Gröbner basis algorithms
   - `multivariate`/`univariate` - Polynomial implementations (supports graded polynomial rings)

4. **`resolution`** - Homological algebra
   - Main Resolution struct for computing projective resolutions

5. **`export`** - Data export functionality for spectral sequences

### Key Design Patterns

- **Generic over grading types**: Most structures are parameterized by grading `G: Grading`
- **Arc-wrapped shared ownership**: Coalgebras are typically wrapped in `Arc<>` for shared references
- **Trait-based architecture**: Core functionality defined through `Comodule` and `ComoduleMorphism` traits
- **Parallel processing**: Uses `rayon` for parallel computation
- **Serialization support**: Structures support `serde` for JSON export

### Input/Output
- Input files are text-based and located in `examples/` directory
- Coalgebras and comodules can be parsed from text files
- Results can be exported to JSON format for spectral sequence visualization
- The main binary processes examples from `examples/polynomial/A.txt`

### Test Organization
- Tests are organized in `tests/` subdirectories within each module
- Integration tests are in the root `tests/` directory
- Test modules follow the pattern `<module_name>_tests.rs`

## Development Notes

- The codebase is actively evolving with both k-linear (`kcomodule`) and R-linear (`rcomodule`) implementations
- New files like `tensor.rs`, `groups.rs`, `module.rs`, `ring.rs` suggest ongoing expansion
- Pay attention to grading constraints - many operations require `OrderedGrading`
- Memory efficiency is important for large computations - consider using `FlatMatrix` for dense operations
- The library supports homological algebra over both finite fields and graded polynomial rings over fields
- R-linear morphisms (`rmorphism`) handle computations over polynomial rings, complementing the k-linear versions