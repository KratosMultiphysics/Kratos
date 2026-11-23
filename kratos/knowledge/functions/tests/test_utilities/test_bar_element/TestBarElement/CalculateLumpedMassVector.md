---
type: C++ Method
title: CalculateLumpedMassVector
resource: tests/test_utilities/test_bar_element.cpp#L336-L365
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestBarElement::CalculateLumpedMassVector( VectorType& rLumpedMassVector, const ProcessInfo& rCurrentProcessInfo) const`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [X0](../../../../../functions/includes/node/X0.md)
- [Y0](../../../../../functions/includes/node/Y0.md)