---
type: Python Method
title: __PrepareOutputSettings
resource: python_scripts/orchestrators/sequential_orchestrator.py#L143-L194
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
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/GetNumberOfStages
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __PrepareOutputSettings(self, stage_name: str) -> str:`

# Calls

- [GetProject](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [IsBool](../../../../../functions/sources/kratos_parameters/Parameters/IsBool.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [IsString](../../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [__GetExecutionList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)
- [GetNumberOfStages](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/GetNumberOfStages.md)
- [__GetStagesCheckpointFolder](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder.md)
- [SetString](../../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [AddString](../../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [remove](../../../../../functions/sources/kratos_filesystem/remove.md)
- [SetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/SetStringArray.md)
- [RemoveValue](../../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)