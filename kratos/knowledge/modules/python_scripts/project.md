---
type: Python Module
title: project
resource: python_scripts/project.py#L1-L136
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/sys
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pickle
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
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
  - target: external/kratosmultiphysics-analysis-stage
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Project](../../classes/python_scripts/project/Project.md)
- [__init__](../../functions/python_scripts/project/Project/__init__.md)
- [GetModel](../../functions/python_scripts/project/Project/GetModel.md)
- [GetSettings](../../functions/python_scripts/project/Project/GetSettings.md)
- [GetOutputData](../../functions/python_scripts/project/Project/GetOutputData.md)
- [GetActiveStages](../../functions/python_scripts/project/Project/GetActiveStages.md)
- [AddActiveStage](../../functions/python_scripts/project/Project/AddActiveStage.md)
- [RemoveActiveStage](../../functions/python_scripts/project/Project/RemoveActiveStage.md)
- [Save](../../functions/python_scripts/project/Project/Save.md)
- [Load](../../functions/python_scripts/project/Project/Load.md)

# Imports

- `sys`
- `pickle`
- `pathlib`
- `importlib`
- `typing`
- `KratosMultiphysics`
- `KratosMultiphysics.analysis_stage`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)