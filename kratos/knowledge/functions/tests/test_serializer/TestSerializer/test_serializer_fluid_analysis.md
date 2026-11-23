---
type: Python Method
title: test_serializer_fluid_analysis
resource: tests/test_serializer.py#L60-L68
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_serializer/TestSerializer/_prepare_fluid_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_serializer_fluid_analysis(self):`

# Calls

- [_prepare_fluid_test](../../../../functions/tests/test_serializer/TestSerializer/_prepare_fluid_test.md)
- [RunSolutionLoop](../../../../functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop.md)