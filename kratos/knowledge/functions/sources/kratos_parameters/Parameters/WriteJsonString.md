---
type: C++ Method
title: WriteJsonString
resource: sources/kratos_parameters.cpp#L407-L410
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/IsEmpty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/replace_elements_and_condition_process/ReplaceElementsAndConditionsProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/Append
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/InternalSetValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/specialized_spatial_search_factory/SpecializedSpatialSearchFactory/SpecializedSpatialSearchFactory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_kratos_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_kratos_include_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_kratos_more_levels_include_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_vector_json_with_include_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_iterators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_file_serialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_stream_serialization
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const std::string Parameters::WriteJsonString() const`

# Called by

- [IsEmpty](../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/IsEmpty.md)
- [ReplaceElementsAndConditionsProcess](../../../../functions/processes/replace_elements_and_condition_process/ReplaceElementsAndConditionsProcess.md)
- [Clone](../../../../functions/sources/kratos_parameters/Parameters/Clone.md)
- [Append](../../../../functions/sources/kratos_parameters/Parameters/Append.md)
- [save](../../../../functions/sources/kratos_parameters/Parameters/save.md)
- [InternalSetValue](../../../../functions/sources/kratos_parameters/Parameters/InternalSetValue.md)
- [SpecializedSpatialSearchFactory](../../../../functions/spatial_containers/specialized_spatial_search_factory/SpecializedSpatialSearchFactory/SpecializedSpatialSearchFactory.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_kratos_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_kratos_parameters.md)
- [test_kratos_include_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_kratos_include_parameters.md)
- [test_kratos_more_levels_include_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_kratos_more_levels_include_parameters.md)
- [test_vector_json_with_include_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_vector_json_with_include_parameters.md)
- [test_iterators](../../../../functions/tests/test_kratos_parameters/TestParameters/test_iterators.md)
- [test_file_serialization](../../../../functions/tests/test_kratos_parameters/TestParameters/test_file_serialization.md)
- [test_stream_serialization](../../../../functions/tests/test_kratos_parameters/TestParameters/test_stream_serialization.md)