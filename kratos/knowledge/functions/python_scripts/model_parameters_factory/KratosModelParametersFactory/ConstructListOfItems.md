---
type: Python Method
title: ConstructListOfItems
resource: python_scripts/model_parameters_factory.py#L83-L89
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructItem
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/run_modelers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructListOfItems( self, items_list):`

# Calls

- [ConstructItem](../../../../functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructItem.md)

# Called by

- [_CreateModelers](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers.md)
- [__CreateListOfModelers](../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers.md)
- [__CreateListOfOperations](../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations.md)
- [run_modelers](../../../../functions/tests/test_voxel_modeler/run_modelers.md)