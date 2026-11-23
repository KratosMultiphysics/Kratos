---
type: Python Method
title: __init__
resource: python_scripts/assign_mpcs_to_neighbours_process.py#L9-L19
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, settings): # Check for deprecated settings and update them # The setting 'assign_mpcs_every_time_step' was renamed to 'reform_constraints_at_each_step' # This block ensures backward compatibility by updating the old setting to the new one and issuing a warning.`

# Calls

- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)