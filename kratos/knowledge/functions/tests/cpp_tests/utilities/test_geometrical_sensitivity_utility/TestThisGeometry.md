---
type: C++ Function
title: TestThisGeometry
resource: tests/cpp_tests/utilities/test_geometrical_sensitivity_utility.cpp#L124-L158
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_sensitivity_utility/ShapeParameter/Sequence/Sequence
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/CheckDeterminantOfJacobianSensitivityByFiniteDifference
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/CheckShapeFunctionsGradientSensitivityByFiniteDifference
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestThisGeometry(Geometry<Point>& rGeom, GeometryData::IntegrationMethod ThisMethod, double StepSize = 1e-7, double Tolerance = 1e-7)`

# Calls

- [Sequence](../../../../../functions/utilities/geometrical_sensitivity_utility/ShapeParameter/Sequence/Sequence.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [CheckDeterminantOfJacobianSensitivityByFiniteDifference](../../../../../functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/CheckDeterminantOfJacobianSensitivityByFiniteDifference.md)
- [CheckShapeFunctionsGradientSensitivityByFiniteDifference](../../../../../functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/CheckShapeFunctionsGradientSensitivityByFiniteDifference.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/KRATOS_TEST_CASE_IN_SUITE.md)