---
type: C++ Function
title: SumNonHistoricalNodeScalarVariable
resource: utilities/variable_utils.h#L1233-L1254
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_sum_variable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`[[nodiscard]] double SumNonHistoricalNodeScalarVariable( const TVarType& rVar, const ModelPart& rModelPart )`

# Called by

- [test_Sum](../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Sum.md)
- [test_sum_variable](../../../functions/tests/test_variable_utils/TestVariableUtils/test_sum_variable.md)