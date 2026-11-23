---
type: C++ Function
title: AssembleContainerContributions
resource: utilities/sensitivity_builder.cpp#L166-L214
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleContainerContributions( TContainerType& rContainer, AdjointResponseFunction& rResponseFunction, SensitivityBuilderScheme& rSensitivityBuilderScheme, TProxyType& rProxy, const Variable<TDataType>& rDesignVariable, const ProcessInfo& rProcessInfo)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)