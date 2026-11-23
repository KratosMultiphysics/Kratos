---
type: C++ Function
title: CreateSubModelPartHierarchy
resource: mpi/utilities/distributed_model_part_initializer.cpp#L56-L62
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/SplitStringByDelimiter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateSubModelPartHierarchy(ModelPart& rModelPart, const std::string& rModelPartHierarchy)`

# Calls

- [SplitStringByDelimiter](../../../../functions/utilities/string_utilities/SplitStringByDelimiter.md)
- [RecursiveCreateModelParts](../../../../functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts.md)

# Called by

- [CopySubModelPartStructure](../../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure.md)