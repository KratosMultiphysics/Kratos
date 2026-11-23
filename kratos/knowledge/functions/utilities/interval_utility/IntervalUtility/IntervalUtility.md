---
type: C++ Method
title: IntervalUtility
resource: utilities/interval_utility.cpp#L19-L33
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_flag_process/AssignFlagProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IntervalUtility::IntervalUtility(Parameters Settings)`

# Calls

- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [__init__](../../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_flag_process/AssignFlagProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__.md)