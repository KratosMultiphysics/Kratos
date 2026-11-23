---
type: C++ Method
title: TetrahedralMeshOrientationCheck
resource: processes/tetrahedral_mesh_orientation_check.cpp#L41-L52
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
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

`TetrahedralMeshOrientationCheck::TetrahedralMeshOrientationCheck( ModelPart& rModelPart, const bool ThrowErrors, const Flags Options ): Process(), mrModelPart(rModelPart), mOptions(Options)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)