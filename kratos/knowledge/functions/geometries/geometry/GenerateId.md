---
type: C++ Function
title: GenerateId
resource: geometries/geometry.h#L1066-L1079
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/SetIdGeneratedFromString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/SetIdNotSelfAssigned
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/geometry/Geometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/SetId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mesh/Mesh/pGetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mesh/Mesh/RemoveGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mesh/Mesh/HasGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_geometry_container
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline IndexType GenerateId(const std::string& rName)`

# Calls

- [SetIdGeneratedFromString](../../../functions/geometries/geometry/SetIdGeneratedFromString.md)
- [SetIdNotSelfAssigned](../../../functions/geometries/geometry/SetIdNotSelfAssigned.md)

# Called by

- [Geometry](../../../functions/geometries/geometry/Geometry.md)
- [SetId](../../../functions/geometries/geometry/SetId.md)
- [pGetGeometry](../../../functions/includes/mesh/Mesh/pGetGeometry.md)
- [RemoveGeometry](../../../functions/includes/mesh/Mesh/RemoveGeometry.md)
- [HasGeometry](../../../functions/includes/mesh/Mesh/HasGeometry.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_geometry_container](../../../functions/tests/test_model_part/TestModelPart/test_geometry_container.md)