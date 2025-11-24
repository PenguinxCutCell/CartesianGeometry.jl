# CartesianGeometry

![CI](https://github.com/PenguinxCutCell/CartesianGeometry.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/CartesianGeometry.jl/branch/main/graph/badge.svg)

CartesianGeometry.jl is a lightweight Julia library providing utilities for Cartesian-grid geometry operations. The package collects mesh and geometry helpers used by numerical methods that operate on structured Cartesian meshes, including routines for implicit integration, interfaces/borders, VOF initialization, and mesh utilities.

**Features**
- **Implicit integration:** fast, robust helpers for integrating fields over cut cells and geometric domains (see `src/implicitint.jl`).
- **Mesh helpers:** utilities for constructing and querying Cartesian meshes (`src/mesh.jl`).
- **Border/interface handling:** routines for border and interface geometry (`src/border.jl`).
- **VOF initialization:** helpers for initializing volume-of-fluid representations (`src/vofinit.jl`).
- **Utilities:** misc helpers in `src/utils.jl` used across algorithms.

**Quick Install (from GitHub)**
Clone and develop locally:

```bash
git clone https://github.com/PenguinxCutCell/CartesianGeometry.jl.git
cd CartesianGeometry.jl
# Start Julia with the package environment and run tests / develop
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Or add directly from the repository URL inside a Julia REPL:

```julia
] add https://github.com/PenguinxCutCell/CartesianGeometry.jl.git
```

**Running Tests**
- Run the full test suite from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```
- Or run the package test driver directly:

```bash
julia --project=. test/runtests.jl
```

**Usage (overview)**
- Load the package in Julia and call the modules/functions you need. The package entrypoint and primary code lives in `src/` and is organized across files like `implicitint.jl`, `mesh.jl`, `border.jl`, `vofinit.jl`, and `utils.jl`.
- Example (REPL):

```julia
julia> import Pkg; Pkg.activate(".") # if working from repo
julia> using CartesianGeometry

# Explore functions provided in the package modules
# (See source files in src/ for available functions and examples.)
```

**Examples: integrate outputs**

- The `integrate` helper can return different groups of geometric quantities depending on the first argument (a `Tuple` of selectors). Common usages:

- Cell-based quantities (volume, barycentre, interface length, cell types):

```julia
V, bary, interface_length, cell_types = integrate(Tuple{0}, levelset, xyz, T, nan)
```

- Face-based quantities (face areas / integrated face properties):

```julia
As = integrate(Tuple{1}, levelset, xyz, T, nan)
```

- You can also request staggered / centroid-based outputs by passing the computed `bary` as an extra argument:

```julia
Ws = integrate(Tuple{0}, levelset, xyz, T, nan, bary)   # staggered volumes (Ws)
Bs = integrate(Tuple{1}, levelset, xyz, T, nan, bary)   # centroid-face areas (Bs)
```

Variables returned:
- `V` : cell volumes
- `bary` : cell barycentres (centroids)
- `interface_length` : interface length / area inside a cell
- `cell_types` : an encoding of cell type (e.g., cut, full, empty)
- `As` : face areas / face-integrated values
- `Ws` : staggered volumes
- `Bs` : centroid face-area-like quantities


**Repository Structure**
- `src/` : package source files (`CartesianGeometry.jl`, `mesh.jl`, `implicitint.jl`, `border.jl`, `vofinit.jl`, `utils.jl`, etc.).
- `test/`: test scripts and `runtests.jl` used by the test harness.
- `bench/`: example/benchmark scripts.
- `Project.toml` / `Manifest.toml`: package environment.
- `LICENSE`: license for the project.


**License**
- See the `LICENSE` file in the repository root for licensing details.

