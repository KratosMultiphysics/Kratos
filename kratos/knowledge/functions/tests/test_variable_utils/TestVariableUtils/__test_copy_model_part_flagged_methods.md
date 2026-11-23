---
type: Python Method
title: __test_copy_model_part_flagged_methods
resource: tests/test_variable_utils.py#L1135-L1202
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_variable_utils/test_case
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_nodal_non_historical_var_to_non_historical_var
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_element_var
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_condition_var
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __test_copy_model_part_flagged_methods( self, method, get_entity_container, origin_entity_value_setter, origin_entity_value_getter, destination_entity_value_setter, destination_entity_value_getter, add_origin_variables_to_solution_step_variables_list = True, add_destination_variables_to_solution_step_variables_list = True): ##set the origin model part`

# Called by

- [test_case](../../../../functions/tests/test_variable_utils/test_case.md)
- [test_copy_model_part_flagged_nodal_non_historical_var_to_non_historical_var](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_nodal_non_historical_var_to_non_historical_var.md)
- [test_copy_model_part_flagged_element_var](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_element_var.md)
- [test_copy_model_part_flagged_condition_var](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_copy_model_part_flagged_condition_var.md)