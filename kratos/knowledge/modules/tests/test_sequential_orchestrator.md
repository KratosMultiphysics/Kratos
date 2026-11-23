---
type: Python Module
title: test_sequential_orchestrator
resource: tests/test_sequential_orchestrator.py#L1-L91
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/importlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-json-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-project
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-python-solver
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

- [EmptySolver](../../classes/tests/test_sequential_orchestrator/EmptySolver.md)
- [__init__](../../functions/tests/test_sequential_orchestrator/EmptySolver/__init__.md)
- [GetDefaultParameters](../../functions/tests/test_sequential_orchestrator/EmptySolver/GetDefaultParameters.md)
- [ImportModelPart](../../functions/tests/test_sequential_orchestrator/EmptySolver/ImportModelPart.md)
- [GetComputingModelPart](../../functions/tests/test_sequential_orchestrator/EmptySolver/GetComputingModelPart.md)
- [EmptyAnalysisStage](../../classes/tests/test_sequential_orchestrator/EmptyAnalysisStage.md)
- [__init__](../../functions/tests/test_sequential_orchestrator/EmptyAnalysisStage/__init__.md)
- [_CreateSolver](../../functions/tests/test_sequential_orchestrator/EmptyAnalysisStage/_CreateSolver.md)
- [KeepAdvancingSolutionLoop](../../functions/tests/test_sequential_orchestrator/EmptyAnalysisStage/KeepAdvancingSolutionLoop.md)
- [EmptyOperation](../../classes/tests/test_sequential_orchestrator/EmptyOperation.md)
- [__init__](../../functions/tests/test_sequential_orchestrator/EmptyOperation/__init__.md)
- [Execute](../../functions/tests/test_sequential_orchestrator/EmptyOperation/Execute.md)
- [TestSequentialOrchestrator](../../classes/tests/test_sequential_orchestrator/TestSequentialOrchestrator.md)
- [test_sequential_orchestrator](../../functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/test_sequential_orchestrator.md)
- [tearDown](../../functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/tearDown.md)

# Imports

- `pathlib`
- `importlib`
- `KratosMultiphysics`
- `KratosMultiphysics.KratosUnittest`
- `KratosMultiphysics.json_utilities`
- `KratosMultiphysics.kratos_utilities`
- `KratosMultiphysics.project`
- `KratosMultiphysics.python_solver`
- `KratosMultiphysics.analysis_stage`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)