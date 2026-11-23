---
type: C++ Method
title: CalculateNodalSolutionStepSensitivities
resource: utilities/sensitivity_builder.cpp#L515-L534
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SensitivityBuilder::CalculateNodalSolutionStepSensitivities( const std::vector<std::string>& rVariables, ModelPart& rModelPart, AdjointResponseFunction& rResponseFunction, double ScalingFactor)`

# Calls

- [GetVariableLists](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists.md)