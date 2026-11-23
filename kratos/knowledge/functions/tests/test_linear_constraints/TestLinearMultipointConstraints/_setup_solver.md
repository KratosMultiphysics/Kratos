---
type: Python Method
title: _setup_solver
resource: tests/test_linear_constraints.py#L71-L108
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_basic_setup_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_setup_test
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _setup_solver(self, solving_with = "Block", linear_solver = "AMGCL"): #define a minimal newton raphson solver`

# Called by

- [_basic_setup_test](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_basic_setup_test.md)
- [_advanced_setup_test](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_setup_test.md)