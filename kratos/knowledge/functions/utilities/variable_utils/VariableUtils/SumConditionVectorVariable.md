---
type: C++ Method
title: SumConditionVectorVariable
resource: utilities/variable_utils.cpp#L75-L92
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_sum_variable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double, 3> VariableUtils::SumConditionVectorVariable( const ArrayVarType& rVar, const ModelPart& rModelPart )`

# Called by

- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [test_sum_variable](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_sum_variable.md)