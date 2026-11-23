---
type: C++ Method
title: CopySubModelPartStructure
resource: mpi/utilities/distributed_model_part_initializer.cpp#L66-L89
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/utilities/distributed_model_part_initializer/GetSubModelPartHierarchy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DistributedModelPartInitializer::CopySubModelPartStructure()`

# Calls

- [GetSubModelPartHierarchy](../../../../../functions/mpi/utilities/distributed_model_part_initializer/GetSubModelPartHierarchy.md)
- [Broadcast](../../../../../functions/includes/data_communicator/Broadcast.md)
- [CreateSubModelPartHierarchy](../../../../../functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy.md)