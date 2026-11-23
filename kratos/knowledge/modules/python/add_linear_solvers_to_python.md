---
type: C++ Module
title: add_linear_solvers_to_python
resource: python/add_linear_solvers_to_python.cpp#L1-L264
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-add-linear-solvers-to-python-h
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
  - target: external/includes-dof-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-complex-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-linear-solver-ublas-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-reorderer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-direct-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-skyline-lu-factorization-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-skyline-lu-custom-scalar-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-scaling-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-fallback-linear-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-preconditioner-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-diagonal-preconditioner-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-ilu0-preconditioner-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-ilu-preconditioner-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-power-iteration-eigenvalue-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-power-iteration-highest-eigenvalue-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-rayleigh-quotient-iteration-eigenvalue-solver-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AddLinearSolversToPython](../../functions/python/add_linear_solvers_to_python/AddLinearSolversToPython.md)

# Imports

- `"includes/define_python.h"`
- `"python/add_linear_solvers_to_python.h"`
- `"linear_solvers/cg_solver.h"`
- `"linear_solvers/deflated_cg_solver.h"`
- `"linear_solvers/bicgstab_solver.h"`
- `"linear_solvers/tfqmr_solver.h"`
- `"includes/dof.h"`
- `"includes/ublas_complex_interface.h"`
- `"linear_solvers/linear_solver_ublas.h"`
- `"linear_solvers/reorderer.h"`
- `"linear_solvers/direct_solver.h"`
- `"linear_solvers/skyline_lu_factorization_solver.h"`
- `"linear_solvers/skyline_lu_custom_scalar_solver.h"`
- `"linear_solvers/scaling_solver.h"`
- `"linear_solvers/fallback_linear_solver.h"`
- `"linear_solvers/preconditioner.h"`
- `"linear_solvers/diagonal_preconditioner.h"`
- `"linear_solvers/ilu0_preconditioner.h"`
- `"linear_solvers/ilu_preconditioner.h"`
- `"linear_solvers/power_iteration_eigenvalue_solver.h"`
- `"linear_solvers/power_iteration_highest_eigenvalue_solver.h"`
- `"linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)