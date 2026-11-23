---
type: C++ Method
title: Exception
resource: sources/exception.cpp#L22-L27
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/exception/Exception/update_what
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/add_to_call_stack
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Exception::Exception() : std::exception(), mMessage("Unknown Error") , mCallStack()`

# Calls

- [update_what](../../../../functions/sources/exception/Exception/update_what.md)
- [add_to_call_stack](../../../../functions/sources/exception/Exception/add_to_call_stack.md)