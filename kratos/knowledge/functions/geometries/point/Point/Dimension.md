---
type: C++ Method
title: Dimension
resource: geometries/point.h#L175-L178
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/mesh/Mesh/WorkingSpaceDimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateDivisionNumbers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_point/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static constexpr IndexType Dimension()`

# Called by

- [WorkingSpaceDimension](../../../../functions/includes/mesh/Mesh/WorkingSpaceDimension.md)
- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)
- [CalculateBoundingBox](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundingBox.md)
- [CalculateDivisionNumbers](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateDivisionNumbers.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_point/KRATOS_TEST_CASE_IN_SUITE.md)