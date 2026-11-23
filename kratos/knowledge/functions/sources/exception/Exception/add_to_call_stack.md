---
type: C++ Method
title: add_to_call_stack
resource: sources/exception.cpp#L59-L63
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/exception/Exception/update_what
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/exception/Exception/Exception
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Exception::add_to_call_stack(CodeLocation const& TheLocation)`

# Calls

- [update_what](../../../../functions/sources/exception/Exception/update_what.md)

# Called by

- [Exception](../../../../functions/sources/exception/Exception/Exception.md)