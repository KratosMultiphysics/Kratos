---
type: Python Method
title: __init__
resource: python_scripts/assign_vector_by_direction_to_entity_process.py#L23-L170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
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
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/VariableUtils/SumConditionVectorVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [CalculateNormals](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals.md)
- [SumConditionVectorVariable](../../../../functions/utilities/variable_utils/VariableUtils/SumConditionVectorVariable.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [IsNumber](../../../../functions/sources/kratos_parameters/Parameters/IsNumber.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)