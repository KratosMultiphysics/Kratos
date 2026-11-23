---
type: Python Method
title: _test_bdf_time_discretization
resource: tests/test_time_discretization.py#L22-L36
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_generic_BDF
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF4
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF5
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF6
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _test_bdf_time_discretization(self, bdf, exp_results, *delta_times):`

# Calls

- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)

# Called by

- [test_generic_BDF](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_generic_BDF.md)
- [test_BDF1](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF1.md)
- [test_BDF2](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF2.md)
- [test_BDF3](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF3.md)
- [test_BDF4](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF4.md)
- [test_BDF5](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF5.md)
- [test_BDF6](../../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_BDF6.md)