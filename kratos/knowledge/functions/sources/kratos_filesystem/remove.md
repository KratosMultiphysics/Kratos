---
type: C++ Function
title: remove
resource: sources/kratos_filesystem.cpp#L57-L60
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/KratosUnittest/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/WorkFolderScope/__exit__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_post_results_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/cell/Cell/Remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_key_plane_generation_with_refinement/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_levelset_convection/TestLevelSetConvection/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_print_info_in_file/TestPrintNodalInfoInFile/test_nodal_print
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_print_info_in_file/TestPrintElementalInfoInFile/test_elemental_values
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool remove(const std::string& rPath)`

# Called by

- [main](../../../functions/python_scripts/KratosUnittest/main.md)
- [__exit__](../../../functions/python_scripts/KratosUnittest/WorkFolderScope/__exit__.md)
- [Generate](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/Generate.md)
- [_generateClassFrom](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [_generateVariablesFrom](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateVariablesFrom.md)
- [__remove_list_files](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files.md)
- [__remove_post_results_files](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_post_results_files.md)
- [DeleteFileIfExisting](../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)
- [__PrepareOutputSettings](../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)
- [main](../../../functions/python_scripts/testing/run_cpp_mpi_tests/main.md)
- [main](../../../functions/python_scripts/testing/run_cpp_tests/main.md)
- [WriteProperties](../../../functions/sources/model_part_io/ModelPartIO/WriteProperties.md)
- [Remove](../../../functions/spatial_containers/cell/Cell/Remove.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/modeler/test_key_plane_generation_with_refinement/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_model_part_io/KRATOS_TEST_CASE_IN_SUITE.md)
- [tearDown](../../../functions/tests/test_levelset_convection/TestLevelSetConvection/tearDown.md)
- [test_nodal_print](../../../functions/tests/test_print_info_in_file/TestPrintNodalInfoInFile/test_nodal_print.md)
- [test_elemental_values](../../../functions/tests/test_print_info_in_file/TestPrintElementalInfoInFile/test_elemental_values.md)