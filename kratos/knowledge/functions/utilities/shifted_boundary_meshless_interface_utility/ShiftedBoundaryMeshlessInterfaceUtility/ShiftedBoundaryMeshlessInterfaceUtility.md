---
type: C++ Method
title: ShiftedBoundaryMeshlessInterfaceUtility
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L126-L182
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ShiftedBoundaryMeshlessInterfaceUtility::ShiftedBoundaryMeshlessInterfaceUtility( Model& rModel, Parameters ThisParameters)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)