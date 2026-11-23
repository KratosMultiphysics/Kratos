---
type: C++ Module
title: p_multigrid_builder_and_solver
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L1-L1296
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-p-multigrid-builder-and-solver-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-constraint-assembler-factory-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-sparse-utilities-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-diagonal-scaling-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-p-grid-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-p-multigrid-status-stream-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-linear-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-components-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-proxies-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-profiler-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/optional
    resolved_by: tree-sitter
    confidence: exact
  - target: external/variant
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [LinearSystem](../../../../classes/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/LinearSystem.md)
- [Impl](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/Impl.md)
- [ExecuteMultigridLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop.md)
- [ExecuteConstraintLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop.md)
- [Solve](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/Solve.md)
- [CollectDoFs](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/CollectDoFs.md)
- [MakeLhsTopology](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/MakeLhsTopology.md)
- [Assemble](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/Assemble.md)
- [constexpr](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/constexpr.md)
- [PMultigridBuilderAndSolver](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/PMultigridBuilderAndSolver.md)
- [PMultigridBuilderAndSolver](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/PMultigridBuilderAndSolver-2.md)
- [Create](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Create.md)
- [FindNodeIndex](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/FindNodeIndex.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet.md)
- [SetUpSystem](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpSystem.md)
- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ResizeAndInitializeVectors.md)
- [InitializeSolutionStep](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/InitializeSolutionStep.md)
- [FinalizeSolutionStep](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse, TDense>/FinalizeSolutionStep.md)
- [Build](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Build.md)
- [BuildLHS](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildLHS.md)
- [BuildRHS](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildRHS.md)
- [BuildRHSAndSolve](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildRHSAndSolve.md)
- [ApplyConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyConstraints.md)
- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyDirichletConditions.md)
- [BuildAndSolve](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildAndSolve.md)
- [SystemSolve](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SystemSolve.md)
- [CalculateReactions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/CalculateReactions.md)
- [AssignSettings](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/AssignSettings.md)
- [Clear](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Clear.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/GetDefaultParameters.md)
- [Info](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Info.md)
- [GetEquationSystemSize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/GetEquationSystemSize.md)
- [ProjectGrid](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid.md)

# Imports

- `"solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler_factory.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/p_grid.hpp"`
- `"solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp"`
- `"includes/model_part.h"`
- `"spaces/ublas_space.h"`
- `"linear_solvers/linear_solver.h"`
- `"factories/linear_solver_factory.h"`
- `"includes/kratos_components.h"`
- `"utilities/proxies.h"`
- `"utilities/profiler.h"`
- `"utilities/builtin_timer.h"`
- `<optional>`
- `<variant>`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)