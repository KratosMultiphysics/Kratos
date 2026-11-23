---
type: Python Method
title: __GetStagesCheckpointFolder
resource: python_scripts/orchestrators/sequential_orchestrator.py#L133-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetStagesCheckpointFolder(self) -> str:`

# Calls

- [GetProject](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)
- [__PrepareOutputSettings](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)