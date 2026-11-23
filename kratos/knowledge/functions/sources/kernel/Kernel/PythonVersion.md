---
type: C++ Method
title: PythonVersion
resource: sources/kernel.cpp#L202-L204
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/logger_output/LoggerOutput/WriteHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/PrintInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string Kernel::PythonVersion()`

# Called by

- [WriteHeader](../../../../functions/input_output/logger_output/LoggerOutput/WriteHeader.md)
- [PrintInfo](../../../../functions/sources/kernel/Kernel/PrintInfo.md)