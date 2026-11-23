---
type: C++ Function
title: SetFlag
resource: utilities/variable_utils.h#L915-L929
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Element>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Condition>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<MasterSlaveConstraint>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__RemoveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_flag_process/AssignFlagProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__RemoveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__InitializeAssignEntityDerivativesToNodesValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_flag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/InitializeNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetFlag( const Flags& rFlag, const bool FlagValue, TContainerType& rContainer )`

# Called by

- [ResetModelPart](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart.md)
- [ResetModelPart](../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart.md)
- [RemoveModelPartEntities](../../../functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities.md)
- [AddIntersectedElementsModelPartElements](../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements.md)
- [InitializeGradientRecoveryModelPart](../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [Execute](../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute.md)
- [Execute](../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Element>/Execute.md)
- [Execute](../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Condition>/Execute.md)
- [Execute](../../../functions/processes/entity_erase_process/EntitiesEraseProcess<MasterSlaveConstraint>/Execute.md)
- [FillAuxiliaryModelPart](../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart.md)
- [Execute](../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)
- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [ExecuteInitializeSolutionStep](../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/ExecuteInitializeSolutionStep.md)
- [__RemoveConstraints](../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__RemoveConstraints.md)
- [ExecuteInitializeSolutionStep](../../../functions/python_scripts/assign_flag_process/AssignFlagProcess/ExecuteInitializeSolutionStep.md)
- [__RemoveConstraints](../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__RemoveConstraints.md)
- [SetupModelPart](../../../functions/tests/test_coordinate_transformation_utils/SetupModelPart.md)
- [setUp](../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUp.md)
- [__InitializeAssignEntityDerivativesToNodesValues](../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__InitializeAssignEntityDerivativesToNodesValues.md)
- [test_set_flag](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_flag.md)
- [RemoveElementAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings.md)
- [RemoveElementsAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings.md)
- [RemoveConditionAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings.md)
- [RemoveConditionsAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings.md)
- [RemoveOrphanNodesFromSubModelParts](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts.md)
- [AddNeighbours](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours.md)
- [InitializeNormals](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/InitializeNormals.md)
- [ComputeUnitNormalsFromAreaNormals](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals.md)
- [AssignEntityDerivativesToNodes](../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)