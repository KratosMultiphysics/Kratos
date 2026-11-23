---
type: Python Method
title: setUp
resource: tests/test_tensor_adaptors.py#L56-L90
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/ClearNonHistoricalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/value_setter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/flag_setter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/fixity_setter
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def setUp(self):`

# Calls

- [ClearNonHistoricalData](../../../../functions/utilities/variable_utils/ClearNonHistoricalData.md)
- [value_setter](../../../../functions/tests/test_tensor_adaptors/value_setter.md)
- [flag_setter](../../../../functions/tests/test_tensor_adaptors/flag_setter.md)
- [fixity_setter](../../../../functions/tests/test_tensor_adaptors/fixity_setter.md)