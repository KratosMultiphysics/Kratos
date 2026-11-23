---
type: C++ Function
title: SavePointer
resource: includes/serializer.h#L1105-L1125
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/IsDerived
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/serializer/save
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SavePointer(std::string const & rTag, const TDataType * pValue)`

# Calls

- [IsDerived](../../../functions/includes/serializer/IsDerived.md)

# Called by

- [save](../../../functions/includes/serializer/save.md)