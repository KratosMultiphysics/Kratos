---
type: C++ Function
title: FlipFlag
resource: utilities/variable_utils.h#L957-L969
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/Flip
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_flag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FlipFlag( const Flags& rFlag, TContainerType& rContainer )`

# Calls

- [Flip](../../../functions/containers/flags/Flip.md)

# Called by

- [test_set_flag](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_flag.md)