---
type: Python Method
title: test_with_real_model
resource: mpi/tests/test_distributed_model_part_initializer.py#L12-L53
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/DistributedModelPartInitializer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_with_real_model(self):`

# Calls

- [DistributedModelPartInitializer](../../../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/DistributedModelPartInitializer.md)
- [GlobalNumberOfNodes](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [NumberOfSubModelParts](../../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)