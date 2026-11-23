---
type: C++ Method
title: EdgeBasedGradientRecoveryProcess
resource: processes/edge_based_gradient_recovery_process.h#L232-L240
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`EdgeBasedGradientRecoveryProcess( Model& rModel, typename TLinearSolver::Pointer pLinearSolver, Parameters ThisParameters) : EdgeBasedGradientRecoveryProcess(rModel, ThisParameters)`

# Calls

- [InitializeGradientRecoveryStrategy](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryStrategy.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)