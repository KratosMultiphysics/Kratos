---
type: C++ Method
title: DeregisterMappers
resource: sources/kratos_application.cpp#L433-L454
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/registry_item/RegistryItem/KeyConstBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/KeyConstEnd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/kratos_application/KratosApplication/DeregisterApplication
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KratosApplication::DeregisterMappers()`

# Calls

- [KeyConstBegin](../../../../functions/sources/registry_item/RegistryItem/KeyConstBegin.md)
- [KeyConstEnd](../../../../functions/sources/registry_item/RegistryItem/KeyConstEnd.md)

# Called by

- [DeregisterApplication](../../../../functions/sources/kratos_application/KratosApplication/DeregisterApplication.md)