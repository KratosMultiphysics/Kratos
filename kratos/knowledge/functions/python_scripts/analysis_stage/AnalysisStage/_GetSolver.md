---
type: Python Method
title: _GetSolver
resource: python_scripts/analysis_stage.py#L258-L262
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Finalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/GetComputingModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/FinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AdvanceTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AddVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_PrepareModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AddDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_InitializeInternals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_analysis_stage/TestAnalysisStage/testNullSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_analysis_stage/TestAnalysisStage/testRealSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_analysis_stage/TestAnalysisStage/testSolverIsCached
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetSolver(self):`

# Called by

- [Finalize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Finalize.md)
- [GetComputingModelPart](../../../../functions/python_scripts/analysis_stage/AnalysisStage/GetComputingModelPart.md)
- [InitializeSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep.md)
- [SolveSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/SolveSolutionStep.md)
- [FinalizeSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/FinalizeSolutionStep.md)
- [Check](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Check.md)
- [Clear](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Clear.md)
- [_AdvanceTime](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AdvanceTime.md)
- [_AddVariables](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AddVariables.md)
- [_PrepareModelPart](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_PrepareModelPart.md)
- [_AddDofs](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AddDofs.md)
- [_InitializeInternals](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_InitializeInternals.md)
- [testNullSolver](../../../../functions/tests/test_analysis_stage/TestAnalysisStage/testNullSolver.md)
- [testRealSolver](../../../../functions/tests/test_analysis_stage/TestAnalysisStage/testRealSolver.md)
- [testSolverIsCached](../../../../functions/tests/test_analysis_stage/TestAnalysisStage/testSolverIsCached.md)