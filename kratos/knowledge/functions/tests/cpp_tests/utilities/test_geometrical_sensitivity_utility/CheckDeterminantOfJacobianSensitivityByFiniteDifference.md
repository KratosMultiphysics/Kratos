---
type: C++ Function
title: CheckDeterminantOfJacobianSensitivityByFiniteDifference
resource: tests/cpp_tests/utilities/test_geometrical_sensitivity_utility.cpp#L77-L91
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/TestThisGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckDeterminantOfJacobianSensitivityByFiniteDifference(double DetJ_Deriv, unsigned IntegrationPoint, ShapeParameter Deriv, Geometry<Point>& rGeom, GeometryData::IntegrationMethod ThisMethod, double StepSize = 1e-7, double Tolerance = 1e-7)`

# Called by

- [TestThisGeometry](../../../../../functions/tests/cpp_tests/utilities/test_geometrical_sensitivity_utility/TestThisGeometry.md)