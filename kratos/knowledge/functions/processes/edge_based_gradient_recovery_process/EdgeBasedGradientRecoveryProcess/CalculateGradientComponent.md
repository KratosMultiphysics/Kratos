---
type: C++ Method
title: CalculateGradientComponent
resource: processes/edge_based_gradient_recovery_process.h#L604-L631
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateGradientComponent( const OriginGetFunctionType& rOriginValueGetter, const GradientSetFunctionType& rGradientValueSetter, const std::size_t Component)`

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitializeSolutionStep.md)