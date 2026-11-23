---
type: C++ Method
title: IsInt
resource: sources/kratos_parameters.cpp#L552-L555
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_is_methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsInt() const`

# Called by

- [ApplySettings](../../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings.md)
- [StructuredMeshGeneratorProcess](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess.md)
- [ValidateTheDefaultParameters](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_is_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_is_methods.md)