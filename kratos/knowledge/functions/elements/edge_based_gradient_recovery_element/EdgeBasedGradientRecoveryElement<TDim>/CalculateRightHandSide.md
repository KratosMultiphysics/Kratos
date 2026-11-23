---
type: C++ Method
title: CalculateRightHandSide
resource: elements/edge_based_gradient_recovery_element.cpp#L123-L156
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EdgeBasedGradientRecoveryElement<TDim>::CalculateRightHandSide( VectorType &rRigthHandSideVector, const ProcessInfo &rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)