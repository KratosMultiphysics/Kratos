---
type: Python Function
title: __CheckRegistryListIsInModule
resource: python_interface/python_registry_utilities.py#L30-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_interface/python_registry_utilities/RegisterModelersList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/RegisterOperationsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/RegisterProcessesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/RegisterStagesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/python_registry_utilities/RegisterOrchestratorsList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CheckRegistryListIsInModule(PythonRegistryListModule, ListKeyword):`

# Called by

- [RegisterModelersList](../../../functions/python_interface/python_registry_utilities/RegisterModelersList.md)
- [RegisterOperationsList](../../../functions/python_interface/python_registry_utilities/RegisterOperationsList.md)
- [RegisterProcessesList](../../../functions/python_interface/python_registry_utilities/RegisterProcessesList.md)
- [RegisterStagesList](../../../functions/python_interface/python_registry_utilities/RegisterStagesList.md)
- [RegisterOrchestratorsList](../../../functions/python_interface/python_registry_utilities/RegisterOrchestratorsList.md)