---
type: C++ Function
title: Assemble
resource: containers/distributed_system_vector.h#L322-L345
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/IsLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/LocalId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Assemble( const TVectorType& rVectorInput, const TIndexVectorType& EquationId )`

# Calls

- [IsLocal](../../../functions/containers/distributed_numbering/DistributedNumbering/IsLocal.md)
- [LocalId](../../../functions/containers/distributed_numbering/DistributedNumbering/LocalId.md)