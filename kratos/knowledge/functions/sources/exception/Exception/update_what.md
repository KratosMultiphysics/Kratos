---
type: C++ Method
title: update_what
resource: sources/exception.cpp#L65-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/exception/Exception/Exception
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/append_message
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/add_to_call_stack
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Exception::update_what()`

# Called by

- [Exception](../../../../functions/sources/exception/Exception/Exception.md)
- [append_message](../../../../functions/sources/exception/Exception/append_message.md)
- [add_to_call_stack](../../../../functions/sources/exception/Exception/add_to_call_stack.md)