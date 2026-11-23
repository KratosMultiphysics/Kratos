---
type: C++ Method
title: CalculateRightHandSide
resource: future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities.h#L72-L87
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
---

# Signature

`void CalculateRightHandSide( VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo) override`

# Calls

- [GetGeometry](../../../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../../../../functions/geometries/geometry/PointsNumber.md)