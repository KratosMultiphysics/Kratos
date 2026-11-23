---
type: Python Method
title: __GetExecutionList
resource: python_scripts/orchestrators/sequential_orchestrator.py#L78-L97
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/GetNumberOfStages
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetExecutionList(self) -> list:`

# Calls

- [GetProject](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)
- [GetNumberOfStages](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/GetNumberOfStages.md)
- [__GetStagesToCheckpointList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [__PrepareOutputSettings](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)