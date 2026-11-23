---
type: C++ Method
title: GetVariableLists
resource: utilities/sensitivity_builder.cpp#L615-L636
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNodalSolutionStepSensitivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SensitivityBuilder::TSensitivityVariables SensitivityBuilder::GetVariableLists( const std::vector<std::string>& rVariableNames)`

# Called by

- [SensitivityBuilder](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)
- [CalculateNodalSolutionStepSensitivities](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNodalSolutionStepSensitivities.md)
- [CalculateNonHistoricalSensitivities](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities.md)