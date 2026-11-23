---
type: C++ Function
title: CopyFromContiguousDataArray
resource: utilities/container_io_utils.h#L164-L197
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/JoinValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CopyFromContiguousDataArray( TContainerType& rContainer, const TSpanType& rDataSpan, TIntegerType const * pShapeBegin, TIntegerType const * pShapeEnd, const TReferenceGetterType& rReferenceGetter)`

# Calls

- [JoinValues](../../../functions/utilities/string_utilities/JoinValues.md)