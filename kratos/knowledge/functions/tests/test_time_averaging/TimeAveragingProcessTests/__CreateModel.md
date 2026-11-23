---
type: Python Method
title: __CreateModel
resource: tests/test_time_averaging.py#L543-L589
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CreateModel(self):`

# Calls

- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)

# Called by

- [testTimeAveragingNodalHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalAvgProcess.md)
- [testTimeAveragingNodalNonHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalAvgProcess.md)
- [testTimeAveragingElementalNonHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalAvgProcess.md)
- [testTimeAveragingNodalHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalRMSProcess.md)
- [testTimeAveragingNodalNonHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalRMSProcess.md)
- [testTimeAveragingElementalNonHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalRMSProcess.md)