---
type: C++ Method
title: AddEmptyValue
resource: sources/kratos_parameters.cpp#L471-L477
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/csv_points_output_process/Factory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/line_output_process/LineOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_2d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofArrayUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/_SetUpTestModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_value
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_vector_interface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_string_array_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/Check
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Parameters Parameters::AddEmptyValue(const std::string& rEntry)`

# Called by

- [ImportModelPart](../../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [GetDefaultParameters](../../../../functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/GetDefaultParameters.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__init__.md)
- [Factory](../../../../functions/python_scripts/csv_points_output_process/Factory.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [__init__](../../../../functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/line_output_process/LineOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_naca_0012_calculate_distance_to_skin_2d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_2d.md)
- [test_naca_0012_calculate_distance_to_skin_3d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_3d.md)
- [testDofArrayUtilities](../../../../functions/tests/test_dofs/TestDofs/testDofArrayUtilities.md)
- [_SetUpTestModelPart](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/_SetUpTestModelPart.md)
- [test_embedded_skin_mapping](../../../../functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping.md)
- [Check](../../../../functions/tests/test_ensight_output_process/Check.md)
- [test_add_value](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_value.md)
- [test_vector_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_vector_interface.md)
- [test_matrix_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface.md)
- [test_set_string_array_valid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_string_array_valid.md)
- [setUp](../../../../functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp.md)
- [setUp](../../../../functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp.md)
- [Check](../../../../functions/tests/test_vtk_output_process/Check.md)