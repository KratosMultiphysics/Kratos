---
type: Python Function
title: RegisterStagesList
resource: python_interface/python_registry_utilities.py#L22-L24
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/python_registry_utilities/__CheckRegistryListIsInModule
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/__RegisterItemList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/python_registry_utilities/RegisterAll
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RegisterStagesList(PythonModuleName, PythonRegistryListModule):`

# Calls

- [__CheckRegistryListIsInModule](../../../functions/python_interface/python_registry_utilities/__CheckRegistryListIsInModule.md)
- [__RegisterItemList](../../../functions/python_interface/python_registry_utilities/__RegisterItemList.md)

# Called by

- [RegisterAll](../../../functions/python_interface/python_registry_utilities/RegisterAll.md)