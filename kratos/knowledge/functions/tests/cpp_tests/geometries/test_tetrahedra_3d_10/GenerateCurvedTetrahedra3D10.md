---
type: C++ Function
title: GenerateCurvedTetrahedra3D10
resource: tests/cpp_tests/geometries/test_tetrahedra_3d_10.cpp#L54-L67
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateReferenceNodes3D10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Geometry<NodeType>::Pointer GenerateCurvedTetrahedra3D10()`

# Calls

- [GenerateReferenceNodes3D10](../../../../../functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/GenerateReferenceNodes3D10.md)
- [X0](../../../../../functions/includes/node/X0.md)
- [Y0](../../../../../functions/includes/node/Y0.md)
- [Z0](../../../../../functions/includes/node/Z0.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_tetrahedra_3d_10/KRATOS_TEST_CASE_IN_SUITE.md)