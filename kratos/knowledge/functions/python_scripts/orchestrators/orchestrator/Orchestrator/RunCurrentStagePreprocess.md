---
type: Python Method
title: RunCurrentStagePreprocess
resource: python_scripts/orchestrators/orchestrator.py#L122-L146
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunCurrentStagePreprocess(self, stage_name: str, data: Optional[Dict] = None):`

# Calls

- [__CreateListOfModelers](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers.md)
- [__CreateListOfOperations](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)