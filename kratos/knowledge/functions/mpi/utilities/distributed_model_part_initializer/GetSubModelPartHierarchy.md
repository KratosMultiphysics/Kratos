---
type: C++ Function
title: GetSubModelPartHierarchy
resource: mpi/utilities/distributed_model_part_initializer.cpp#L33-L42
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetSubModelPartHierarchy(const ModelPart& rModelPart, std::string& rModelPartHierarchy)`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [CopySubModelPartStructure](../../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure.md)