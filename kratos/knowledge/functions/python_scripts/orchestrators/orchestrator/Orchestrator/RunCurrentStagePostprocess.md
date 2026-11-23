---
type: Python Method
title: RunCurrentStagePostprocess
resource: python_scripts/orchestrators/orchestrator.py#L148-L162
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunCurrentStagePostprocess(self, stage_name: str, data: Optional[Dict] = None):`

# Calls

- [__CreateListOfOperations](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations.md)

# Called by

- [Run](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/Run.md)