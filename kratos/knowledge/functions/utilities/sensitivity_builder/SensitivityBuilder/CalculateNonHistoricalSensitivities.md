---
type: C++ Method
title: CalculateNonHistoricalSensitivities
resource: utilities/sensitivity_builder.cpp#L553-L567
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

`void SensitivityBuilder::CalculateNonHistoricalSensitivities( const std::vector<std::string>& rVariables, ModelPart::ElementsContainerType& rElements, AdjointResponseFunction& rResponseFunction, const ProcessInfo& rProcessInfo, double ScalingFactor)`

# Calls

- [GetVariableLists](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists.md)