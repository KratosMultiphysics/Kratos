---
type: Python Method
title: _advanced_setup_test
resource: tests/test_linear_constraints.py#L333-L363
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_setup_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_add_dofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_set_and_fill_buffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_apply_mpc_constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_setup_solver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_check_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_reset
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_MPC_Constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_LM_MPC_Constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_Double_LM_MPC_Constraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _advanced_setup_test(self, solving_with = "Block", linear_solver = "AMGCL"):`

# Calls

- [_advanced_setup_model_part](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_setup_model_part.md)
- [_add_dofs](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_add_dofs.md)
- [_set_and_fill_buffer](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_set_and_fill_buffer.md)
- [_advanced_apply_mpc_constraints](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_apply_mpc_constraints.md)
- [_setup_solver](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_setup_solver.md)
- [_solve](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_solve.md)
- [_advanced_check_results](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_advanced_check_results.md)
- [_reset](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/_reset.md)

# Called by

- [test_advanced_MPC_Constraints](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_MPC_Constraints.md)
- [test_advanced_LM_MPC_Constraints](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_LM_MPC_Constraints.md)
- [test_advanced_Double_LM_MPC_Constraints](../../../../functions/tests/test_linear_constraints/TestLinearMultipointConstraints/test_advanced_Double_LM_MPC_Constraints.md)