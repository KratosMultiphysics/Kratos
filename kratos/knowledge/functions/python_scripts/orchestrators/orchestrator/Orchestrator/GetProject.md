---
type: Python Method
title: GetProject
resource: python_scripts/orchestrators/orchestrator.py#L164-L167
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetProject(self) -> Project:`

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)
- [__GetExecutionList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)
- [__GetStagesToCheckpointList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [__GetStagesCheckpointFolder](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder.md)
- [__PrepareOutputSettings](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)