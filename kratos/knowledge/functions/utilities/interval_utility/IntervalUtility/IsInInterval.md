---
type: C++ Method
title: IsInInterval
resource: utilities/interval_utility.cpp#L51-L55
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_flag_process/AssignFlagProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/ExecuteFinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/IsOutputStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/SetInitialStateProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IntervalUtility::IsInInterval(double Time)`

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/ExecuteInitializeSolutionStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_flag_process/AssignFlagProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteFinalizeSolutionStep](../../../../functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/ExecuteFinalizeSolutionStep.md)
- [IsOutputStep](../../../../functions/python_scripts/point_output_process/PointOutputProcess/IsOutputStep.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/set_initial_state_process/SetInitialStateProcess/ExecuteInitializeSolutionStep.md)