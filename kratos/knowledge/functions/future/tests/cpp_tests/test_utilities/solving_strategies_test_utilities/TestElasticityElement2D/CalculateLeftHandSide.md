---
type: C++ Method
title: CalculateLeftHandSide
resource: future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities.h#L146-L204
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
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

`void CalculateLeftHandSide( MatrixType &rLeftHandSideMatrix, const ProcessInfo &rCurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../../../../functions/geometries/geometry/PointsNumber.md)
- [X0](../../../../../../../functions/includes/node/X0.md)
- [Y0](../../../../../../../functions/includes/node/Y0.md)