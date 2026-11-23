---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/processes/test_structured_mesh_generator_process.cpp#L27-L78
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(StructuredMeshGeneratorProcessHexahedra, KratosCoreFastSuite)`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)