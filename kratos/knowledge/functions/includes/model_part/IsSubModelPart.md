---
type: C++ Function
title: IsSubModelPart
resource: includes/model_part.h#L1873-L1876
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/InsertEntityRange
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/InsertEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CopyCommonData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateCommunicators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/CopyCommonData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/AddEntitiesFromIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/Reset
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CloneSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CloneTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/OverwriteSolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/ReduceTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AssignNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodeFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/SetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddTable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveTableFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/pGetProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemovePropertiesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElementFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddMasterSlaveConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewMasterSlaveConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveGeometryFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/SetBufferSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddElementsWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddConditionsWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongingsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AddElementsWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AddConditionsWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliaryAddEntitiesWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_utils/AddNodesFromOrderedContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: |-
      functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
          TEntityType, TContainerType>/BooleanOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsSubModelPart() const`

# Called by

- [FullName](../../../functions/includes/model_part/FullName.md)
- [InsertEntityRange](../../../functions/includes/model_part/InsertEntityRange.md)
- [InsertEntities](../../../functions/includes/model_part/InsertEntities.md)
- [CopyCommonData](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CopyCommonData.md)
- [CreateCommunicators](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CreateCommunicators.md)
- [PopulateCommunicators](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators.md)
- [CopyCommonData](../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/CopyCommonData.md)
- [AddEntitiesFromIds](../../../functions/sources/model_part/AddEntitiesFromIds.md)
- [Reset](../../../functions/sources/model_part/ModelPart/Reset.md)
- [CloneSolutionStep](../../../functions/sources/model_part/ModelPart/CloneSolutionStep.md)
- [CloneTimeStep](../../../functions/sources/model_part/ModelPart/CloneTimeStep.md)
- [CreateTimeStep](../../../functions/sources/model_part/ModelPart/CreateTimeStep.md)
- [OverwriteSolutionStepData](../../../functions/sources/model_part/ModelPart/OverwriteSolutionStepData.md)
- [ReduceTimeStep](../../../functions/sources/model_part/ModelPart/ReduceTimeStep.md)
- [AddNode](../../../functions/sources/model_part/ModelPart/AddNode.md)
- [AssignNode](../../../functions/sources/model_part/ModelPart/AssignNode.md)
- [RemoveNodeFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveNodeFromAllLevels.md)
- [SetNodalSolutionStepVariablesList](../../../functions/sources/model_part/ModelPart/SetNodalSolutionStepVariablesList.md)
- [AddTable](../../../functions/sources/model_part/ModelPart/AddTable.md)
- [RemoveTableFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveTableFromAllLevels.md)
- [AddProperties](../../../functions/sources/model_part/ModelPart/AddProperties.md)
- [RecursivelyHasProperties](../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [pGetProperties](../../../functions/sources/model_part/ModelPart/pGetProperties.md)
- [RemovePropertiesFromAllLevels](../../../functions/sources/model_part/ModelPart/RemovePropertiesFromAllLevels.md)
- [AddElement](../../../functions/sources/model_part/ModelPart/AddElement.md)
- [RemoveElementFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveElementFromAllLevels.md)
- [AddMasterSlaveConstraint](../../../functions/sources/model_part/ModelPart/AddMasterSlaveConstraint.md)
- [CreateNewMasterSlaveConstraint](../../../functions/sources/model_part/ModelPart/CreateNewMasterSlaveConstraint.md)
- [RemoveMasterSlaveConstraintFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveMasterSlaveConstraintFromAllLevels.md)
- [AddCondition](../../../functions/sources/model_part/ModelPart/AddCondition.md)
- [RemoveConditionFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveConditionFromAllLevels.md)
- [AddGeometry](../../../functions/sources/model_part/ModelPart/AddGeometry.md)
- [CreateNewGeometry](../../../functions/sources/model_part/ModelPart/CreateNewGeometry.md)
- [RemoveGeometryFromAllLevels](../../../functions/sources/model_part/ModelPart/RemoveGeometryFromAllLevels.md)
- [SetBufferSize](../../../functions/sources/model_part/ModelPart/SetBufferSize.md)
- [PrintData](../../../functions/sources/model_part/ModelPart/PrintData.md)
- [save](../../../functions/sources/model_part/ModelPart/save.md)
- [load](../../../functions/sources/model_part/ModelPart/load.md)
- [AddElementsWithNodes](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddElementsWithNodes.md)
- [AddConditionsWithNodes](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddConditionsWithNodes.md)
- [RemoveElementAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongingsFromAllLevels.md)
- [RemoveConditionAndBelongingsFromAllLevels](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongingsFromAllLevels.md)
- [AddElementsWithNodes](../../../functions/utilities/auxiliar_model_part_utilities/AddElementsWithNodes.md)
- [AddConditionsWithNodes](../../../functions/utilities/auxiliar_model_part_utilities/AddConditionsWithNodes.md)
- [AuxiliaryAddEntitiesWithNodes](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliaryAddEntitiesWithNodes.md)
- [AddNodesFromOrderedContainer](../../../functions/utilities/model_part_utils/AddNodesFromOrderedContainer.md)
- [BooleanOperation](../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>/BooleanOperation.md)