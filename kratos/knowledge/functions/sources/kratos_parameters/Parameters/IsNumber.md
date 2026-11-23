---
type: C++ Method
title: IsNumber
resource: sources/kratos_parameters.cpp#L536-L539
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/components_to_functions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling/Scaling/Scaling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsNumber() const`

# Called by

- [ApplyConstantVectorValueProcess](../../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__.md)
- [components_to_functions](../../../../functions/python_scripts/set_initial_state_process/components_to_functions.md)
- [Scaling](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling/Scaling/Scaling.md)