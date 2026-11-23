---
type: C++ Module
title: implicit_scheme
resource: future/solving_strategies/schemes/implicit_scheme.h#L1-L1635
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/execution
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-dof-utilities-dof-array-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-entities-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-openmp-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/future-solving-strategies-builders-builder-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/future-solving-strategies-builders-block-builder-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/future-solving-strategies-builders-elimination-builder-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ImplicitThreadLocalStorage](../../../../classes/future/solving_strategies/schemes/implicit_scheme/ImplicitThreadLocalStorage.md)
- [ImplicitScheme](../../../../classes/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme.md)
- [ImplicitScheme](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ImplicitScheme.md)
- [ImplicitScheme](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ImplicitScheme-2.md)
- [ImplicitScheme](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ImplicitScheme-3.md)
- [Create](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Create.md)
- [Clone](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Clone.md)
- [Initialize](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Initialize.md)
- [InitializeSolutionStep](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeSolutionStep.md)
- [Predict](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Predict.md)
- [FinalizeSolutionStep](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/FinalizeSolutionStep.md)
- [InitializeNonLinIteration](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeNonLinIteration.md)
- [FinalizeNonLinIteration](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/FinalizeNonLinIteration.md)
- [SetUpDofArrays](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpDofArrays.md)
- [SetUpSystemIds](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpSystemIds.md)
- [Build](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build.md)
- [BuildWithSafeAssemble](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble.md)
- [BuildWithThreadLocal](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal.md)
- [BuildWithLocalAllocation](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation.md)
- [Build](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build-2.md)
- [Build](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build-3.md)
- [BuildMassMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMassMatrix.md)
- [BuildDampingMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildDampingMatrix.md)
- [BuildMasterSlaveConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints.md)
- [BuildLinearSystemConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildLinearSystemConstraints.md)
- [ApplyLinearSystemConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ApplyLinearSystemConstraints.md)
- [CalculateReactions](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateReactions.md)
- [Update](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Update.md)
- [UpdateConstraintsOnlyDofs](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/UpdateConstraintsOnlyDofs.md)
- [CalculateUpdateVector](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateUpdateVector.md)
- [Clear](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Clear.md)
- [Check](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Check.md)
- [GetDefaultParameters](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetDefaultParameters.md)
- [Name](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Name.md)
- [SetMoveMesh](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetMoveMesh.md)
- [SetReformDofsAtEachStep](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetReformDofsAtEachStep.md)
- [SetEchoLevel](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetEchoLevel.md)
- [GetMoveMesh](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetMoveMesh.md)
- [GetReformDofsAtEachStep](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetReformDofsAtEachStep.md)
- [GetEchoLevel](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/GetEchoLevel.md)
- [Info](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Info.md)
- [PrintInfo](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/PrintInfo.md)
- [PrintData](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/PrintData.md)
- [InitializeLinearSystem](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem.md)
- [CalculateLocalSystemContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateLocalSystemContribution.md)
- [CalculateRightHandSideContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateRightHandSideContribution.md)
- [CalculateLeftHandSideContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateLeftHandSideContribution.md)
- [CalculateMassMatrixContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateMassMatrixContribution.md)
- [CalculateDampingMatrixContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateDampingMatrixContribution.md)
- [CalculateConstraintContribution](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateConstraintContribution.md)
- [MoveMesh](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/MoveMesh.md)
- [ValidateAndAssignParameters](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ValidateAndAssignParameters.md)
- [AssignSettings](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/AssignSettings.md)

# Imports

- `<algorithm>`
- `<execution>`
- `"includes/kratos_parameters.h"`
- `"includes/model_part.h"`
- `"utilities/builtin_timer.h"`
- `"utilities/dof_utilities/dof_array_utilities.h"`
- `"utilities/entities_utilities.h"`
- `"utilities/openmp_utils.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/timer.h"`
- `"future/solving_strategies/builders/builder.h"`
- `"future/solving_strategies/builders/block_builder.h"`
- `"future/solving_strategies/builders/elimination_builder.h"`

# Member of

- [KratosMultiphysics](../../../../packages/KratosMultiphysics.md)