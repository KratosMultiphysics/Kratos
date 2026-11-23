---
type: C++ Module
title: edge_based_gradient_recovery_process
resource: processes/edge_based_gradient_recovery_process.h#L1-L681
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-global-nodal-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-residualbased-block-builder-and-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-residualbased-incrementalupdate-static-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-strategies-residualbased-linear-strategy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GradientDataHandler](../../classes/processes/edge_based_gradient_recovery_process/GradientDataHandler.md)
- [GetOrigin](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<double, true>/GetOrigin.md)
- [GetOrigin](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<double, false>/GetOrigin.md)
- [GetOrigin](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<array_1d<double,3>, true>/GetOrigin.md)
- [GetOrigin](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<array_1d<double,3>, false>/GetOrigin.md)
- [SetGradient](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<double, true>/SetGradient.md)
- [SetGradient](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<double, false>/SetGradient.md)
- [SetGradient](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<array_1d<double,3>, true>/SetGradient.md)
- [SetGradient](../../functions/processes/edge_based_gradient_recovery_process/GradientDataHandler<array_1d<double,3>, false>/SetGradient.md)
- [EdgeBasedGradientRecoveryProcess](../../classes/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess.md)
- [EdgeBasedGradientRecoveryProcess](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess.md)
- [EdgeBasedGradientRecoveryProcess](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess-2.md)
- [Check](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/Check.md)
- [Execute](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/Execute.md)
- [ExecuteInitialize](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteFinalizeSolutionStep](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteFinalizeSolutionStep.md)
- [ExecuteFinalize](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteFinalize.md)
- [Clear](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/Clear.md)
- [GetDefaultParameters](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/GetDefaultParameters.md)
- [Info](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/Info.md)
- [PrintInfo](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/PrintInfo.md)
- [PrintData](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/PrintData.md)
- [EdgeBasedGradientRecoveryProcess](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess-3.md)
- [InitializeGradientRecoveryModelPart](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [CheckAndAssignSettings](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CheckAndAssignSettings.md)
- [InitializeGradientRecoveryStrategy](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryStrategy.md)
- [CalculateGradientComponent](../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CalculateGradientComponent.md)

# Imports

- `"includes/define.h"`
- `"processes/find_global_nodal_neighbours_process.h"`
- `"solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"`
- `"solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"`
- `"solving_strategies/strategies/residualbased_linear_strategy.h"`
- `"utilities/variable_utils.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)