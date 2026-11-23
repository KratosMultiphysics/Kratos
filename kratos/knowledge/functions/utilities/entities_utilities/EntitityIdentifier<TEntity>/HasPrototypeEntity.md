---
type: C++ Method
title: HasPrototypeEntity
resource: utilities/entities_utilities.cpp#L66-L69
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_entities_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool EntitityIdentifier<TEntity>::HasPrototypeEntity(const GeometryType& rGeometry) const`

# Called by

- [CreateEntitiesFromGeometries](../../../../functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometries.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_entities_utilities/KRATOS_TEST_CASE_IN_SUITE.md)