---
type: C++ Method
title: DeleteModelPart
resource: containers/model.cpp#L64-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_linear_master_slave_constraint/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Model::DeleteModelPart( const std::string& rModelPartName )`

# Calls

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)

# Called by

- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)
- [GenerateIntersectedEdgesElementsModelPart](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)
- [Execute](../../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [InitializeGradientRecoveryModelPart](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [Clear](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Clear.md)
- [ReGenerateDistanceModelPart](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart.md)
- [SetupModelPart](../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_linear_master_slave_constraint/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestIODebug](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)