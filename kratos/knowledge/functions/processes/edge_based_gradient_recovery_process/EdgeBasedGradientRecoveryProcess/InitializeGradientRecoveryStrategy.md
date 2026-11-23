---
type: C++ Method
title: InitializeGradientRecoveryStrategy
resource: processes/edge_based_gradient_recovery_process.h#L582-L602
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeGradientRecoveryStrategy(BuilderAndSolverPointerType pBuilderAndSolver)`

# Calls

- [InitializeGradientRecoveryModelPart](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [EdgeBasedGradientRecoveryProcess](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess.md)