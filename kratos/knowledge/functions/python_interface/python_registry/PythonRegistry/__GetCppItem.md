---
type: Python Method
title: __GetCppItem
resource: python_interface/python_registry.py#L268-L281
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/python_registry/PythonRegistry/__GetCppRegistryGetFunctionNameFromItemKeyword
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/python_registry/PythonRegistry/__getitem__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry/PythonRegistry/keys
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry/PythonRegistry/values
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetCppItem(self, Name): # Check if current item has value`

# Calls

- [__GetCppRegistryGetFunctionNameFromItemKeyword](../../../../functions/python_interface/python_registry/PythonRegistry/__GetCppRegistryGetFunctionNameFromItemKeyword.md)

# Called by

- [__getitem__](../../../../functions/python_interface/python_registry/PythonRegistry/__getitem__.md)
- [keys](../../../../functions/python_interface/python_registry/PythonRegistry/keys.md)
- [values](../../../../functions/python_interface/python_registry/PythonRegistry/values.md)