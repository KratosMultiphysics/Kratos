---
type: C++ Method
title: CalculateLocalSystem
resource: elements/edge_based_gradient_recovery_element.cpp#L28-L68
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

`void EdgeBasedGradientRecoveryElement<TDim>::CalculateLocalSystem( MatrixType &rLeftHandSideMatrix, VectorType &rRigthHandSideVector, const ProcessInfo &rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)