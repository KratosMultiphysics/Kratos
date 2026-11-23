---
type: C++ Function
title: CheckSameModelPartUsingSkinDistanceProcess
resource: processes/check_same_modelpart_using_skin_distance_process.h#L71-L85
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/check_same_model_part_using_skin_process/Factory
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit CheckSameModelPartUsingSkinDistanceProcess( Model& rModel, Parameters ThisParameters = Parameters(R"({})") ) : mrSkinModelPart1(rModel.GetModelPart(ThisParameters["skin_model_part_1_name"].GetString())), mrSkinModelPart2(rModel.GetModelPart(ThisParameters["skin_model_part_2_name"].GetString())), mThisParameters(ThisParameters)`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyValidateAndAssignDefaults](../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)

# Called by

- [Factory](../../../functions/python_scripts/check_same_model_part_using_skin_process/Factory.md)