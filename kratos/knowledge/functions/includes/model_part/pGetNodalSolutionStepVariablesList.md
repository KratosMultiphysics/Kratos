---
type: C++ Function
title: pGetNodalSolutionStepVariablesList
resource: includes/model_part.h#L544-L547
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CopyCommonData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/CopyCommonData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`VariablesList::Pointer pGetNodalSolutionStepVariablesList() const`

# Called by

- [CopyCommonData](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CopyCommonData.md)
- [CopyCommonData](../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/CopyCommonData.md)
- [GenerateNode](../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateNode.md)
- [AddEdgeNodes](../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes.md)