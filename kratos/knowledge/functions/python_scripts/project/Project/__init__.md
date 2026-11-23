---
type: Python Method
title: __init__
resource: python_scripts/project.py#L23-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kernel/Kernel/Version
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kernel/Kernel/BuildType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, settings: KratosMultiphysics.Parameters) -> None:`

# Calls

- [Version](../../../../functions/sources/kernel/Kernel/Version.md)
- [BuildType](../../../../functions/sources/kernel/Kernel/BuildType.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)