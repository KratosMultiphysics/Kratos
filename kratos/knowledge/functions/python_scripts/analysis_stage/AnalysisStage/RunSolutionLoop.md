---
type: Python Method
title: RunSolutionLoop
resource: python_scripts/analysis_stage.py#L59-L68
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AdvanceTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/OutputSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/Run
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_serializer/TestSerializer/test_serializer_fluid_analysis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunSolutionLoop(self):`

# Calls

- [_AdvanceTime](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AdvanceTime.md)
- [OutputSolutionStep](../../../../functions/python_scripts/analysis_stage/AnalysisStage/OutputSolutionStep.md)

# Called by

- [Run](../../../../functions/python_scripts/analysis_stage/AnalysisStage/Run.md)
- [test_serializer_fluid_analysis](../../../../functions/tests/test_serializer/TestSerializer/test_serializer_fluid_analysis.md)