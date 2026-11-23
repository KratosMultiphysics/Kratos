---
type: C++ Method
title: EquationIdVector
resource: elements/edge_based_gradient_recovery_element.cpp#L159-L173
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EdgeBasedGradientRecoveryElement<2>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDofPosition](../../../../functions/includes/node/GetDofPosition.md)