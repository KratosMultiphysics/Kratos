---
type: Python Method
title: __GetStagesToCheckpointList
resource: python_scripts/orchestrators/sequential_orchestrator.py#L99-L131
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __GetStagesToCheckpointList(self) -> list:`

# Calls

- [GetProject](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [IsBool](../../../../../functions/sources/kratos_parameters/Parameters/IsBool.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [__GetExecutionList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)
- [IsStringArray](../../../../../functions/sources/kratos_parameters/Parameters/IsStringArray.md)
- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)