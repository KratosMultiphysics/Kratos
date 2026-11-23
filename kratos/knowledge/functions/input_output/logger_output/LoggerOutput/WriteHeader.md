---
type: C++ Method
title: WriteHeader
resource: input_output/logger_output.cpp#L42-L51
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/Version
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/OSName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/PythonVersion
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/Compiler
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LoggerOutput::WriteHeader()`

# Calls

- [Version](../../../../functions/sources/kernel/Kernel/Version.md)
- [OSName](../../../../functions/sources/kernel/Kernel/OSName.md)
- [PythonVersion](../../../../functions/sources/kernel/Kernel/PythonVersion.md)
- [Compiler](../../../../functions/sources/kernel/Kernel/Compiler.md)