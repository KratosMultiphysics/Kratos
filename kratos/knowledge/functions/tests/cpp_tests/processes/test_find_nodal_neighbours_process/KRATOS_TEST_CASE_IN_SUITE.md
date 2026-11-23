---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/processes/test_find_nodal_neighbours_process.cpp#L33-L60
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/CreateSphereTriangularMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(FindNodalNeighboursProcess1, KratosCoreFastSuite)`

# Calls

- [Create2DGeometry](../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry.md)
- [Create3DGeometry](../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry.md)
- [CreateSphereTriangularMesh](../../../../../functions/tests/test_utilities/cpp_tests_utilities/CreateSphereTriangularMesh.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)