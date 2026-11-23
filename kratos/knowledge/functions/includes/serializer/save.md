---
type: C++ Function
title: save
resource: includes/serializer.h#L598-L603
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/save_trace_point
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/save_associative_container
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/IsDerived
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/SavePointer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void save(std::string const & rTag, std::array<TDataType, TDataSize> const& rObject)`

# Calls

- [save_trace_point](../../../functions/includes/serializer/save_trace_point.md)
- [save_associative_container](../../../functions/includes/serializer/save_associative_container.md)
- [IsDerived](../../../functions/includes/serializer/IsDerived.md)
- [SavePointer](../../../functions/includes/serializer/SavePointer.md)