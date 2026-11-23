---
type: C++ Method
title: RemoveNodes
resource: sources/model_part.cpp#L499-L532
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/RemoveCommonNodesFromSlaveModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels
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
---

# Signature

`void ModelPart::RemoveNodes(Flags IdentifierFlag)`

# Called by

- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)
- [CollapseNodes](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes.md)
- [RemoveCommonNodesFromSlaveModelPart](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/RemoveCommonNodesFromSlaveModelPart.md)
- [Execute](../../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute.md)
- [Execute](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/Execute.md)
- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveElementAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings.md)
- [RemoveElementsAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings.md)
- [RemoveConditionAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings.md)
- [RemoveConditionsAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings.md)
- [RemoveOrphanNodesFromSubModelParts](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts.md)