---
type: C++ Function
title: AssignFlags
resource: containers/flags.h#L178-L182
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/ReduceValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/Clear
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignFlags(Flags const& rOther)`

# Called by

- [ReduceValues](../../../functions/mpi/includes/mpi_communicator/ReduceValues.md)
- [Clear](../../../functions/sources/model_part/ModelPart/Clear.md)