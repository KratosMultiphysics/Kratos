---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_geometry_utils.cpp#L99-L113
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTriangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTetrahedra
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/SideLenghts2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTriangle3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTriangleDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/JacobianOnInitialConfiguration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/ShapeFunctionsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/DeformationGradient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SetSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/EvaluateHistoricalVariableValueAtGaussPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsPointDistanceToTriangleOutOfPlane, KratosCoreFastSuite)`

# Calls

- [GenerateExampleTriangle](../../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTriangle.md)
- [GenerateExampleTetrahedra](../../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTetrahedra.md)
- [GenerateExampleLine](../../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleLine.md)
- [SideLenghts2D](../../../../../functions/utilities/geometry_utilities/SideLenghts2D.md)
- [GenerateExampleTriangle3D](../../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/GenerateExampleTriangle3D.md)
- [CalculateTriangleDistances](../../../../../functions/utilities/geometry_utilities/CalculateTriangleDistances.md)
- [GetDefaultIntegrationMethod](../../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [JacobianOnInitialConfiguration](../../../../../functions/utilities/geometry_utilities/JacobianOnInitialConfiguration.md)
- [ShapeFunctionsGradients](../../../../../functions/utilities/geometry_utilities/ShapeFunctionsGradients.md)
- [DeformationGradient](../../../../../functions/utilities/geometry_utilities/DeformationGradient.md)
- [SetSolutionStepVariablesList](../../../../../functions/includes/node/SetSolutionStepVariablesList.md)
- [EvaluateHistoricalVariableValueAtGaussPoint](../../../../../functions/utilities/geometry_utilities/GeometryUtils/EvaluateHistoricalVariableValueAtGaussPoint.md)
- [ProjectedIsInside](../../../../../functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside.md)