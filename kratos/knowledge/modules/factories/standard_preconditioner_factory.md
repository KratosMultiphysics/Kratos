---
type: C++ Module
title: standard_preconditioner_factory
resource: factories/standard_preconditioner_factory.cpp#L1-L53
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-standard-preconditioner-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
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
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RegisterPreconditioners](../../functions/factories/standard_preconditioner_factory/RegisterPreconditioners.md)

# Imports

- `"includes/define.h"`
- `"factories/standard_preconditioner_factory.h"`
- `"spaces/ublas_space.h"`
- `"linear_solvers/preconditioner.h"`
- `"linear_solvers/diagonal_preconditioner.h"`
- `"linear_solvers/ilu0_preconditioner.h"`
- `"linear_solvers/ilu_preconditioner.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)