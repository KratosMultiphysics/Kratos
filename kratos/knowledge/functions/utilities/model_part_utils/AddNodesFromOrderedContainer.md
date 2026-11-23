---
type: C++ Function
title: AddNodesFromOrderedContainer
resource: utilities/model_part_utils.h#L110-L130
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_utils/JoinOrderedNodesContainerType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/AddNodesToModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void AddNodesFromOrderedContainer(ModelPart& rTargetModelPart, TIteratorType iNodesBegin, TIteratorType iNodesEnd)`

# Calls

- [JoinOrderedNodesContainerType](../../../functions/utilities/model_part_utils/JoinOrderedNodesContainerType.md)
- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [AddNodesToModelPart](../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/AddNodesToModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE.md)