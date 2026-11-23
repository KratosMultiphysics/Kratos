---
type: Python Method
title: test_set_flag
resource: tests/test_variable_utils.py#L384-L465
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/FlipFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ResetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_set_flag(self):`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [FlipFlag](../../../../functions/utilities/variable_utils/FlipFlag.md)
- [ResetFlag](../../../../functions/utilities/variable_utils/ResetFlag.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)