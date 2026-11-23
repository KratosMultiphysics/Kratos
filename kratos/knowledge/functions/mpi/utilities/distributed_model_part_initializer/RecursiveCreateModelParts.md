---
type: C++ Function
title: RecursiveCreateModelParts
resource: mpi/utilities/distributed_model_part_initializer.cpp#L44-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/SplitStringByDelimiter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RecursiveCreateModelParts(ModelPart& rModelPart, const std::string& rModelPartName)`

# Calls

- [SplitStringByDelimiter](../../../../functions/utilities/string_utilities/SplitStringByDelimiter.md)
- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [CreateSubModelPartHierarchy](../../../../functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy.md)