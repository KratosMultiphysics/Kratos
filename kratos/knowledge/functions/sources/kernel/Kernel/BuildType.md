---
type: C++ Method
title: BuildType
resource: sources/kernel.cpp#L189-L191
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_version/GetBuildType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/project/Project/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string Kernel::BuildType()`

# Calls

- [GetBuildType](../../../../functions/sources/kratos_version/GetBuildType.md)

# Called by

- [__init__](../../../../functions/python_scripts/project/Project/__init__.md)