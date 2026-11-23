---
type: Python Function
title: Compute_RHS_and_LHS
resource: python_scripts/sympy_fe_utilities.py#L498-L510
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/Compute_RHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/Compute_LHS
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications=False):`

# Calls

- [Compute_RHS](../../../functions/python_scripts/sympy_fe_utilities/Compute_RHS.md)
- [Compute_LHS](../../../functions/python_scripts/sympy_fe_utilities/Compute_LHS.md)