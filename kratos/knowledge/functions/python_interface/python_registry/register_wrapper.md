---
type: Python Function
title: register_wrapper
resource: python_interface/python_registry.py#L322-L357
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry/PythonRegistry/AddItem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def register_wrapper(Class): # Get the list of compiled applications # This will be used for checking the available submodules`

# Calls

- [GetListOfAvailableApplications](../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [AddItem](../../../functions/python_interface/python_registry/PythonRegistry/AddItem.md)