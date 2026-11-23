---
type: C++ Function
title: save_associative_container
resource: includes/serializer.h#L1145-L1154
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/save_trace_point
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/serializer/save
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void save_associative_container(std::string const & rTag, TMapType const& rObject)`

# Calls

- [save_trace_point](../../../functions/includes/serializer/save_trace_point.md)

# Called by

- [save](../../../functions/includes/serializer/save.md)