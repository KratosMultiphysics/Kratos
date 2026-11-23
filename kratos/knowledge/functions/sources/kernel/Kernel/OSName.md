---
type: C++ Method
title: OSName
resource: sources/kernel.cpp#L198-L200
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_version/GetOSName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/logger_output/LoggerOutput/WriteHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/PrintInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string Kernel::OSName()`

# Calls

- [GetOSName](../../../../functions/sources/kratos_version/GetOSName.md)

# Called by

- [WriteHeader](../../../../functions/input_output/logger_output/LoggerOutput/WriteHeader.md)
- [PrintInfo](../../../../functions/sources/kernel/Kernel/PrintInfo.md)