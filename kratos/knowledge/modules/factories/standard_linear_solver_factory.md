---
type: C++ Module
title: standard_linear_solver_factory
resource: factories/standard_linear_solver_factory.cpp#L1-L109
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/factories-standard-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-cg-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-deflated-cg-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-bicgstab-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-tfqmr-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-amgcl-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-amgcl-ns-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-scaling-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-fallback-linear-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-monotonicity-preserving-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-skyline-lu-custom-scalar-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RegisterLinearSolvers](../../functions/factories/standard_linear_solver_factory/RegisterLinearSolvers.md)
- [RegisterLinearSolvers](../../functions/factories/standard_linear_solver_factory/RegisterLinearSolvers-2.md)

# Imports

- `"factories/standard_linear_solver_factory.h"`
- `"linear_solvers/cg_solver.h"`
- `"linear_solvers/deflated_cg_solver.h"`
- `"linear_solvers/bicgstab_solver.h"`
- `"linear_solvers/tfqmr_solver.h"`
- `"linear_solvers/amgcl_solver.h"`
- `"linear_solvers/amgcl_ns_solver.h"`
- `"linear_solvers/scaling_solver.h"`
- `"linear_solvers/fallback_linear_solver.h"`
- `"linear_solvers/monotonicity_preserving_solver.h"`
- `"linear_solvers/skyline_lu_custom_scalar_solver.h"`
- `"spaces/ublas_space.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)