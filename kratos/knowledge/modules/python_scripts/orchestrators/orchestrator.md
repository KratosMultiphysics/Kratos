---
type: Python Module
title: orchestrator
resource: python_scripts/orchestrators/orchestrator.py#L1-L186
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/abc
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sys
    resolved_by: tree-sitter
    confidence: exact
  - target: external/importlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/typing
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-project
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-analysis-stage
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-model-parameters-factory
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Orchestrator](../../../classes/python_scripts/orchestrators/orchestrator/Orchestrator.md)
- [__init__](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__init__.md)
- [Run](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/Run.md)
- [CheckStageSettings](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/CheckStageSettings.md)
- [CreateStage](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/CreateStage.md)
- [RunCurrentStagePreprocess](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePreprocess.md)
- [RunCurrentStagePostprocess](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePostprocess.md)
- [GetProject](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [__CreateListOfModelers](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers.md)
- [__CreateListOfOperations](../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations.md)

# Imports

- `abc`
- `sys`
- `importlib`
- `typing`
- `KratosMultiphysics`
- `KratosMultiphysics.project`
- `KratosMultiphysics.analysis_stage`
- `KratosMultiphysics.model_parameters_factory`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)