---
type: C++ Method
title: AddEmptyArray
resource: sources/kratos_parameters.cpp#L1101-L1107
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetSolvers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_empty_array
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Write
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddEmptyArray(const std::string& rEntry)`

# Called by

- [GetParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetParameters.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [SetSolvers](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetSolvers.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_empty_array](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_empty_array.md)
- [WriteTagsToJson](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)