---
type: C++ Module
title: builder_and_solver
resource: solving_strategies/builder_and_solvers/builder_and_solver.h#L1-L840
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [BuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuilderAndSolver.md)
- [BuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuilderAndSolver-2.md)
- [BuilderAndSolver](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuilderAndSolver-3.md)
- [Create](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Create.md)
- [GetCalculateReactionsFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetCalculateReactionsFlag.md)
- [SetCalculateReactionsFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetCalculateReactionsFlag.md)
- [GetDofSetIsInitializedFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetDofSetIsInitializedFlag.md)
- [SetDofSetIsInitializedFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetDofSetIsInitializedFlag.md)
- [GetReshapeMatrixFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag.md)
- [SetReshapeMatrixFlag](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetReshapeMatrixFlag.md)
- [GetEquationSystemSize](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetEquationSystemSize.md)
- [GetLinearSystemSolver](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetLinearSystemSolver.md)
- [SetLinearSystemSolver](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetLinearSystemSolver.md)
- [BuildLHS](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildLHS.md)
- [BuildRHS](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildRHS.md)
- [Build](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Build.md)
- [BuildLHS_CompleteOnFreeRows](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildLHS_CompleteOnFreeRows.md)
- [BuildLHS_Complete](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildLHS_Complete.md)
- [SystemSolve](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SystemSolve.md)
- [BuildAndSolve](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildAndSolve.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [BuildRHSAndSolve](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/BuildRHSAndSolve.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ApplyDirichletConditions.md)
- [ApplyDirichletConditions_LHS](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ApplyDirichletConditions_LHS.md)
- [ApplyDirichletConditions_RHS](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ApplyDirichletConditions_RHS.md)
- [ApplyRHSConstraints](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ApplyRHSConstraints.md)
- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ApplyConstraints.md)
- [SetUpDofSet](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetUpDofSet.md)
- [SetUpSystem](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetUpSystem.md)
- [ResizeAndInitializeVectors](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ResizeAndInitializeVectors.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/InitializeSolutionStep.md)
- [FinalizeSolutionStep](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/FinalizeSolutionStep.md)
- [CalculateReactions](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/CalculateReactions.md)
- [Clear](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Clear.md)
- [Check](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Check.md)
- [GetDefaultParameters](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetDefaultParameters.md)
- [Name](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Name.md)
- [SetEchoLevel](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/SetEchoLevel.md)
- [GetEchoLevel](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetEchoLevel.md)
- [Info](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/Info.md)
- [PrintInfo](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/PrintInfo.md)
- [PrintData](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/PrintData.md)
- [ValidateAndAssignParameters](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ValidateAndAssignParameters.md)
- [AssignSettings](../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/AssignSettings.md)

# Imports

- `<set>`
- `"includes/define.h"`
- `"includes/model_part.h"`
- `"solving_strategies/schemes/scheme.h"`
- `"includes/kratos_parameters.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)