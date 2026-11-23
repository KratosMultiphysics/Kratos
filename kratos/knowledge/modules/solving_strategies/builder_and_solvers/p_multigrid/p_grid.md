---
type: C++ Module
title: p_grid
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L1-L658
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-p-grid-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-p-multigrid-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-constraint-assembler-factory-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-sparse-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-status-stream-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-components-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-profiler-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sparse-matrix-multiplication-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/limits
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [PGrid](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/PGrid.md)
- [PGrid](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/PGrid-2.md)
- [PGrid](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/PGrid-3.md)
- [MakeLhsTopology](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/MakeLhsTopology.md)
- [Assemble](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Assemble.md)
- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ApplyDirichletConditions.md)
- [ApplyConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ApplyConstraints.md)
- [Initialize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Initialize.md)
- [ExecuteMultigridLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ExecuteMultigridLoop.md)
- [ExecuteConstraintLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ExecuteConstraintLoop.md)
- [ApplyCoarseCorrection](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ApplyCoarseCorrection.md)
- [Finalize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Finalize.md)
- [Restrict](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Restrict.md)
- [Prolong](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Prolong.md)
- [Clear](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/Clear.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/GetDefaultParameters.md)

# Imports

- `"solving_strategies/builder_and_solvers/p_multigrid/p_grid.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler_factory.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp"`
- `"spaces/ublas_space.h"`
- `"factories/linear_solver_factory.h"`
- `"includes/kratos_components.h"`
- `"utilities/profiler.h"`
- `"utilities/sparse_matrix_multiplication_utility.h"`
- `"utilities/builtin_timer.h"`
- `<limits>`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)