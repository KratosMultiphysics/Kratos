---
type: C++ Method
title: HasModelPart
resource: containers/model.cpp#L178-L202
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/RenameModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/controllers/output_controller/OutputController/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_io_modeler/CadIoModeler/SetupGeometryModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_io_modeler/CadIoModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/ReGenerateConvectionModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/Check
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
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_GetOrCreateSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_connectivity_preserve_modeler/TestConnectivityPreserveModeler/test_setup_with_json_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_delete_model_parts_modeler/TestDeleteModelPartModeler/testDeleteModelPartModelerOneModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_delete_model_parts_modeler/TestDeleteModelPartModeler/testDeleteModelPartModelerMultipleModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_import_obj_modeler/TestImportOBJModeler/test_ImportOBJModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model/TestModel/test_model
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model/TestModel/test_model_serialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model/TestModel/test_model_serialization_with_pickling
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/TestSerializer/test_serializer_loading
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Model::HasModelPart(const std::string& rFullModelPartName) const`

# Calls

- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)
- [RenameModelPart](../../../../functions/containers/model/Model/RenameModelPart.md)
- [Check](../../../../functions/controllers/output_controller/OutputController/Check.md)
- [SetupGeometryModel](../../../../functions/modeler/cad_io_modeler/CadIoModeler/SetupGeometryModel.md)
- [SetupModelPart](../../../../functions/modeler/cad_io_modeler/CadIoModeler/SetupModelPart.md)
- [SetupModelPart](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart.md)
- [SetupModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/SetupModelPart.md)
- [DuplicateMesh](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateMesh.md)
- [SetupModelPart](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/SetupModelPart.md)
- [GenerateIntersectedEdgesElementsModelPart](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)
- [InitializeGradientRecoveryModelPart](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [ReGenerateConvectionModelPart](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/ReGenerateConvectionModelPart.md)
- [Check](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/Check.md)
- [Clear](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Clear.md)
- [ReGenerateDistanceModelPart](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart.md)
- [SetupModelPart](../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/SetupModelPart.md)
- [_GetOrCreateSubModelPart](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_GetOrCreateSubModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_setup_with_json_parameters](../../../../functions/tests/test_connectivity_preserve_modeler/TestConnectivityPreserveModeler/test_setup_with_json_parameters.md)
- [testDeleteModelPartModelerOneModelPart](../../../../functions/tests/test_delete_model_parts_modeler/TestDeleteModelPartModeler/testDeleteModelPartModelerOneModelPart.md)
- [testDeleteModelPartModelerMultipleModelParts](../../../../functions/tests/test_delete_model_parts_modeler/TestDeleteModelPartModeler/testDeleteModelPartModelerMultipleModelParts.md)
- [test_ImportOBJModeler](../../../../functions/tests/test_import_obj_modeler/TestImportOBJModeler/test_ImportOBJModeler.md)
- [test_model](../../../../functions/tests/test_model/TestModel/test_model.md)
- [test_model_serialization](../../../../functions/tests/test_model/TestModel/test_model_serialization.md)
- [test_model_serialization_with_pickling](../../../../functions/tests/test_model/TestModel/test_model_serialization_with_pickling.md)
- [test_serializer_loading](../../../../functions/tests/test_serializer/TestSerializer/test_serializer_loading.md)