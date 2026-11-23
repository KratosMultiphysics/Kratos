---
type: C++ Function
title: ApplyVector
resource: utilities/variable_utils.h#L1193-L1213
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/CheckVariableExists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_apply_vector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyVector( const TVarType& rVar, const Vector& rData, NodesContainerType& rNodes )`

# Calls

- [CheckVariableExists](../../../functions/utilities/variable_utils/CheckVariableExists.md)

# Called by

- [test_apply_vector](../../../functions/tests/test_variable_utils/TestVariableUtils/test_apply_vector.md)