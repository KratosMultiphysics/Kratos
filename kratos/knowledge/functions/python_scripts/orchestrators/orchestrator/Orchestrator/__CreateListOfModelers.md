---
type: Python Method
title: __CreateListOfModelers
resource: python_scripts/orchestrators/orchestrator.py#L169-L175
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePreprocess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CreateListOfModelers(self, stage_name: str) -> list:`

# Calls

- [GetModel](../../../../../functions/python_scripts/project/Project/GetModel.md)
- [ConstructListOfItems](../../../../../functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems.md)

# Called by

- [RunCurrentStagePreprocess](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePreprocess.md)