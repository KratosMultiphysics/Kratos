---
type: C++ Method
title: ExecuteInitialize
resource: processes/edge_based_gradient_recovery_process.h#L286-L302
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExecuteInitialize() override`

# Calls

- [InitializeGradientRecoveryModelPart](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [SetNonHistoricalVariableToZero](../../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)