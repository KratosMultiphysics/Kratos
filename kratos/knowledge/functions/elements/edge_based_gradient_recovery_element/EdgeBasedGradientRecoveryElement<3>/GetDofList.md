---
type: C++ Method
title: GetDofList
resource: elements/edge_based_gradient_recovery_element.cpp#L210-L224
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

`void EdgeBasedGradientRecoveryElement<3>::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)