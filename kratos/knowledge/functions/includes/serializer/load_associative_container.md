---
type: C++ Function
title: load_associative_container
resource: includes/serializer.h#L1130-L1142
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/load_trace_point
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/serializer/load
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void load_associative_container(std::string const & rTag, TMapType& rObject)`

# Calls

- [load_trace_point](../../../functions/includes/serializer/load_trace_point.md)

# Called by

- [load](../../../functions/includes/serializer/load.md)