---
type: C++ Method
title: Version
resource: sources/kernel.cpp#L194-L196
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_version/GetVersionString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/logger_output/LoggerOutput/WriteHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/PrintInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string Kernel::Version()`

# Calls

- [GetVersionString](../../../../functions/sources/kratos_version/GetVersionString.md)

# Called by

- [WriteHeader](../../../../functions/input_output/logger_output/LoggerOutput/WriteHeader.md)
- [__init__](../../../../functions/python_scripts/project/Project/__init__.md)
- [PrintInfo](../../../../functions/sources/kernel/Kernel/PrintInfo.md)