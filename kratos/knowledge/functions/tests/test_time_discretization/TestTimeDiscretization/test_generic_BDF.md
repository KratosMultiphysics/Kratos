---
type: Python Method
title: test_generic_BDF
resource: tests/test_time_discretization.py#L38-L58
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/time_discretization/BDF
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/_get_expected_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/_test_bdf_time_discretization
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_generic_BDF(self):`

# Calls

- [BDF](../../../../functions/utilities/time_discretization/BDF.md)
- [_get_expected_results](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/_get_expected_results.md)
- [_test_bdf_time_discretization](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/_test_bdf_time_discretization.md)