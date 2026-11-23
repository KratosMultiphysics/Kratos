---
type: C++ Function
title: SortedInsert
resource: containers/pointer_vector_set.h#L1130-L1184
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector_set/GetPointer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/pointer_vector_set/insert
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/insert
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SortedInsert(TIteratorType first, TIteratorType last)`

# Calls

- [GetPointer](../../../functions/containers/pointer_vector_set/GetPointer.md)

# Called by

- [insert](../../../functions/containers/pointer_vector_set/insert.md)
- [insert](../../../functions/containers/pointer_vector_set/PointerVectorSet/insert.md)