---
type: C++ Method
title: GetStringArray
resource: sources/kratos_parameters.cpp#L693-L702
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/IntegrationValuesExtrapolationToNodesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/SubModelPartSkinDetectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/TimeAveragingProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/auxiliary_solver_utilities/AddVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/auxiliary_solver_utilities/AddDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/auxiliary_solver_utilities/AddAuxiliaryDofsToDofsWithReactionsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GenerateVariableListFromInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GenerateFlagsListFromInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/CreateVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/AssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/ExecuteBasicEnsightOutputProcessCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_get_string_array_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_get_string_array_invalid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_string_array_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_string_array_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/GetDofsListFromGenericEntitiesSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineFlagsUsed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineTimeIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/CheckCompatibleConstitutiveLaws
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> Parameters::GetStringArray() const`

# Calls

- [IsStringArray](../../../../functions/sources/kratos_parameters/Parameters/IsStringArray.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [ExecuteInitialize](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/ExecuteInitialize.md)
- [IntegrationValuesExtrapolationToNodesProcess](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/IntegrationValuesExtrapolationToNodesProcess.md)
- [SubModelPartSkinDetectionProcess](../../../../functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/SubModelPartSkinDetectionProcess.md)
- [TimeAveragingProcess](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/TimeAveragingProcess.md)
- [__init__](../../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__.md)
- [AddVariables](../../../../functions/python_scripts/auxiliary_solver_utilities/AddVariables.md)
- [AddDofs](../../../../functions/python_scripts/auxiliary_solver_utilities/AddDofs.md)
- [AddAuxiliaryDofsToDofsWithReactionsList](../../../../functions/python_scripts/auxiliary_solver_utilities/AddAuxiliaryDofsToDofsWithReactionsList.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [GenerateVariableListFromInput](../../../../functions/python_scripts/kratos_utilities/GenerateVariableListFromInput.md)
- [GenerateFlagsListFromInput](../../../../functions/python_scripts/kratos_utilities/GenerateFlagsListFromInput.md)
- [SetupModelPart](../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart.md)
- [__GetExecutionList](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)
- [__GetStagesToCheckpointList](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [__init__](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__.md)
- [CreateVariablesList](../../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/CreateVariablesList.md)
- [AssignSettings](../../../../functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/AssignSettings.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_assign_unique_model_part_collection_tag_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [ExecuteBasicEnsightOutputProcessCheck](../../../../functions/tests/test_ensight_output_process/ExecuteBasicEnsightOutputProcessCheck.md)
- [test_get_string_array_valid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_get_string_array_valid.md)
- [test_get_string_array_invalid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_get_string_array_invalid.md)
- [test_set_string_array_valid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_string_array_valid.md)
- [test_add_string_array_valid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_string_array_valid.md)
- [ReadTagsFromJson](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/ReadTagsFromJson.md)
- [SensitivityBuilder](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)
- [GetDofsListFromGenericEntitiesSpecifications](../../../../functions/utilities/specifications_utilities/GetDofsListFromGenericEntitiesSpecifications.md)
- [AddMissingVariablesFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromEntitiesList.md)
- [AddMissingVariablesFromSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications.md)
- [AddMissingDofsFromEntitiesList](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromEntitiesList.md)
- [AddMissingDofsFromSpecifications](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromSpecifications.md)
- [DetermineFlagsUsed](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineFlagsUsed.md)
- [DetermineTimeIntegration](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineTimeIntegration.md)
- [DetermineIfCompatibleGeometries](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries.md)
- [CheckCompatibleConstitutiveLaws](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/CheckCompatibleConstitutiveLaws.md)