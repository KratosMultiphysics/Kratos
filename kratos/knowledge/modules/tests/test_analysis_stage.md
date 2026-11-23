---
type: Python Module
title: test_analysis_stage
resource: tests/test_analysis_stage.py#L1-L95
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics-kratosunittest-as-kratosunittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-analysis-stage
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-python-null-solver
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TestAnalysisStage](../../classes/tests/test_analysis_stage/TestAnalysisStage.md)
- [_AnalysisStageWithDummySolver](../../classes/tests/test_analysis_stage/_AnalysisStageWithDummySolver.md)
- [DummySolver](../../classes/tests/test_analysis_stage/DummySolver.md)
- [__init__](../../functions/tests/test_analysis_stage/DummySolver/__init__.md)
- [AddVariables](../../functions/tests/test_analysis_stage/DummySolver/AddVariables.md)
- [SolveSolutionStep](../../functions/tests/test_analysis_stage/DummySolver/SolveSolutionStep.md)
- [_CreateSolver](../../functions/tests/test_analysis_stage/_AnalysisStageWithDummySolver/_CreateSolver.md)
- [testNullSolver](../../functions/tests/test_analysis_stage/TestAnalysisStage/testNullSolver.md)
- [testRealSolver](../../functions/tests/test_analysis_stage/TestAnalysisStage/testRealSolver.md)
- [testSolverIsCached](../../functions/tests/test_analysis_stage/TestAnalysisStage/testSolverIsCached.md)

# Imports

- `KratosMultiphysics.KratosUnittest as KratosUnittest`
- `KratosMultiphysics`
- `KratosMultiphysics.analysis_stage`
- `KratosMultiphysics.python_null_solver`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)