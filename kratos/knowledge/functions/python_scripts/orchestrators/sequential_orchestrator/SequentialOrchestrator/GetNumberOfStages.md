---
type: Python Method
title: GetNumberOfStages
resource: python_scripts/orchestrators/sequential_orchestrator.py#L73-L76
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetNumberOfStages(self) -> int:`

# Calls

- [__GetExecutionList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)

# Called by

- [__PrepareOutputSettings](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)