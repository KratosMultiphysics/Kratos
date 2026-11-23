---
type: C++ Method
title: IsVector
resource: sources/kratos_parameters.cpp#L602-L613
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataJSON
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_is_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_vector_interface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsVector() const`

# Called by

- [IdentifyDataJSON](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataJSON.md)
- [StructuredMeshGeneratorProcess](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess.md)
- [ValidateTheDefaultParameters](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_is_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_is_methods.md)
- [test_vector_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_vector_interface.md)