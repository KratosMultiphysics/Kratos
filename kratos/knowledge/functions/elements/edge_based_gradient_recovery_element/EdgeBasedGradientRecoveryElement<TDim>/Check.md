---
type: C++ Method
title: Check
resource: elements/edge_based_gradient_recovery_element.cpp#L71-L89
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int EdgeBasedGradientRecoveryElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)