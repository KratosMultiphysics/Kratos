---
type: C++ Module
title: amgcl_solver_impl
resource: future/linear_solvers/amgcl_solver_impl.cpp#L1-L227
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/amgcl-adapter-crs-tuple-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-ublas-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-zero-copy-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-block-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-builtin-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-value-type-static-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-make-solver-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-amg-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-coarsening-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-relaxation-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-solver-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-preconditioner-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-vexcl-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-vexcl-static-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/future-linear-solvers-amgcl-solver-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [register_vexcl_static_matrix_type](../../../functions/future/linear_solvers/amgcl_solver_impl/register_vexcl_static_matrix_type.md)
- [AMGCLScalarSolve](../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLScalarSolve.md)
- [AMGCLBlockSolve](../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLBlockSolve.md)
- [AMGCLSolve](../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLSolve.md)

# Imports

- `<amgcl/adapter/crs_tuple.hpp>`
- `<amgcl/adapter/ublas.hpp>`
- `<amgcl/adapter/zero_copy.hpp>`
- `<amgcl/adapter/block_matrix.hpp>`
- `<amgcl/backend/builtin.hpp>`
- `<amgcl/value_type/static_matrix.hpp>`
- `<amgcl/make_solver.hpp>`
- `<amgcl/amg.hpp>`
- `<amgcl/coarsening/runtime.hpp>`
- `<amgcl/relaxation/runtime.hpp>`
- `<amgcl/solver/runtime.hpp>`
- `<amgcl/preconditioner/runtime.hpp>`
- `<amgcl/backend/vexcl.hpp>`
- `<amgcl/backend/vexcl_static_matrix.hpp>`
- `"future/linear_solvers/amgcl_solver.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)