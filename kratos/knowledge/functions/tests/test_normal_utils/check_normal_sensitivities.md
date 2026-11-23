---
type: Python Function
title: check_normal_sensitivities
resource: tests/test_normal_utils.py#L68-L74
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/get_matrix_row
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_normal_sensitivities(condition, row_index, analytical_normal_shape_sensitivities):`

# Calls

- [assertVectorAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [get_matrix_row](../../../functions/tests/test_normal_utils/get_matrix_row.md)

# Called by

- [FiniteDifferenceNormalShapeSensitivityTest](../../../functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest.md)