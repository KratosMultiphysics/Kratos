---
type: C++ Function
title: SetSolutionStepVariablesList
resource: includes/node.h#L524-L527
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AssignNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/SetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList)`

# Called by

- [GenerateNode](../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateNode.md)
- [AddEdgeNodes](../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes.md)
- [CreateNewNode](../../../functions/sources/model_part/ModelPart/CreateNewNode.md)
- [AssignNode](../../../functions/sources/model_part/ModelPart/AssignNode.md)
- [SetNodalSolutionStepVariablesList](../../../functions/sources/model_part/ModelPart/SetNodalSolutionStepVariablesList.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [AddMissingVariablesFromSpecifications](../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingVariablesFromSpecifications.md)