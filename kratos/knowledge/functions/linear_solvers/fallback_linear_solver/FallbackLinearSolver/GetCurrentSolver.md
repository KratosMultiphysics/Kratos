---
type: C++ Method
title: GetCurrentSolver
resource: linear_solvers/fallback_linear_solver.h#L754-L762
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/AdditionalPhysicalDataIsNeeded
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ProvideAdditionalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/IsConsistent
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/IsNotConsistent
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/UpdateSolverIndex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`LinearSolverPointer GetCurrentSolver()`

# Called by

- [InitializeSolutionStep](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/InitializeSolutionStep.md)
- [FinalizeSolutionStep](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FinalizeSolutionStep.md)
- [PerformSolutionStep](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/PerformSolutionStep.md)
- [Solve](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/Solve.md)
- [AdditionalPhysicalDataIsNeeded](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/AdditionalPhysicalDataIsNeeded.md)
- [ProvideAdditionalData](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ProvideAdditionalData.md)
- [IsConsistent](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/IsConsistent.md)
- [IsNotConsistent](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/IsNotConsistent.md)
- [UpdateSolverIndex](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/UpdateSolverIndex.md)