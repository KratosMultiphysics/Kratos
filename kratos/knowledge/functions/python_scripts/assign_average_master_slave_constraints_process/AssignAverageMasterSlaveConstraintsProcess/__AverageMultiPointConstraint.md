---
type: Python Method
title: __AverageMultiPointConstraint
resource: python_scripts/assign_average_master_slave_constraints_process.py#L62-L84
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __AverageMultiPointConstraint(self, variable):`

# Calls

- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)

# Called by

- [ExecuteInitialize](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitializeSolutionStep.md)