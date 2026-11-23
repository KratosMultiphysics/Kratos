---
type: C++ Method
title: ExecuteInitializeSolutionStep
resource: processes/edge_based_gradient_recovery_process.h#L304-L325
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CalculateGradientComponent
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExecuteInitializeSolutionStep() override`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [CalculateGradientComponent](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CalculateGradientComponent.md)