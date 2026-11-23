---
type: C++ Method
title: SetInt
resource: sources/kratos_parameters.cpp#L834-L837
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/json_io/KratosJsonIO/EnforceConsecutiveOrdering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/line_search_strategy/OverrideDefaultSettingsWithParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
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
  - target: functions/tests/test_kratos_parameters/TestParameters/test_kratos_copy_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::SetInt(const int Value)`

# Called by

- [EnforceConsecutiveOrdering](../../../../functions/includes/json_io/KratosJsonIO/EnforceConsecutiveOrdering.md)
- [Execute](../../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [__init__](../../../../functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [OverrideDefaultSettingsWithParameters](../../../../functions/solving_strategies/strategies/line_search_strategy/OverrideDefaultSettingsWithParameters.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_naca_0012_calculate_distance_to_skin_2d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_2d.md)
- [test_naca_0012_calculate_distance_to_skin_3d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_3d.md)
- [testDofArrayUtilities](../../../../functions/tests/test_dofs/TestDofs/testDofArrayUtilities.md)
- [_SetUpTestModelPart](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/_SetUpTestModelPart.md)
- [test_embedded_skin_mapping](../../../../functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping.md)
- [test_kratos_copy_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_kratos_copy_parameters.md)
- [test_set_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_methods.md)
- [setUp](../../../../functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp.md)
- [setUp](../../../../functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp.md)