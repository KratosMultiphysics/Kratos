---
type: C++ Function
title: Flip
resource: containers/flags.h#L215-L219
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_flags/TestFlags/testFlagFlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/activation_utilities/ActivationUtilities/ActivateElementsAndConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/FlipFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Flip(const Flags ThisFlag)`

# Called by

- [testFlagFlip](../../../functions/tests/test_flags/TestFlags/testFlagFlip.md)
- [ActivateElementsAndConditions](../../../functions/utilities/activation_utilities/ActivationUtilities/ActivateElementsAndConditions.md)
- [FlipFlag](../../../functions/utilities/variable_utils/FlipFlag.md)