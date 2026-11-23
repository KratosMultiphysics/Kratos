---
type: C++ Function
title: JoinValues
resource: utilities/string_utilities.h#L232-L247
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tensor_adaptors/tensor_adaptor_utils/GetTensorShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/container_io_utils/CopyFromContiguousDataArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`[[nodiscard]] std::string JoinValues( TIteratorType Begin, TIteratorType End, const std::string& rDelimiter)`

# Called by

- [GetTensorShape](../../../functions/tensor_adaptors/tensor_adaptor_utils/GetTensorShape.md)
- [CopyFromContiguousDataArray](../../../functions/utilities/container_io_utils/CopyFromContiguousDataArray.md)