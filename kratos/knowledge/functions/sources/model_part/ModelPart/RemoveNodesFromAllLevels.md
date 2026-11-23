---
type: C++ Method
title: RemoveNodesFromAllLevels
resource: sources/model_part.cpp#L534-L538
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/ModelPartRemoveNodesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_remove_nodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/__generate_discontinous_case
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/Clear
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::RemoveNodesFromAllLevels(Flags IdentifierFlag)`

# Calls

- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)

# Called by

- [ResetModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart.md)
- [ResetModelPart](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart.md)
- [test_remove_nodes_parallel_interfaces](../../../../functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces.md)
- [Execute](../../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Node>/Execute.md)
- [ModelPartRemoveNodesFromAllLevels](../../../../functions/python/add_model_part_to_python/ModelPartRemoveNodesFromAllLevels.md)
- [test_remove_nodes](../../../../functions/tests/test_model_part/TestModelPart/test_remove_nodes.md)
- [__generate_discontinous_case](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/__generate_discontinous_case.md)
- [Clear](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/Clear.md)