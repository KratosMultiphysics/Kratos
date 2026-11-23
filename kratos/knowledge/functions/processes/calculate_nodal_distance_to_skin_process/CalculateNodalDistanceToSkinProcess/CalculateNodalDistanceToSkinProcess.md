---
type: C++ Method
title: CalculateNodalDistanceToSkinProcess
resource: processes/calculate_nodal_distance_to_skin_process.cpp#L58-L87
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
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateNodalDistanceToSkinProcess::CalculateNodalDistanceToSkinProcess( Model& rModel, Parameters ThisParameters ) : mrVolumeModelPart(rModel.GetModelPart(ThisParameters["volume_model_part"].GetString())), mrSkinModelPart(rModel.GetModelPart(ThisParameters["skin_model_part"].GetString()))`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)