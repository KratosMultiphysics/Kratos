---
type: C++ Method
title: FromJSONCheckResultProcess
resource: processes/from_json_check_result_process.cpp#L37-L44
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
---

# Signature

`FromJSONCheckResultProcess::FromJSONCheckResultProcess( Model& rModel, Parameters ThisParameters ) : mpModelPart(&[](Model& rModel,const Parameters& rParameters)->ModelPart&{const std::string& r_model_part_name = rParameters["model_part_name"].GetString(); const std::string sub_model_part_name = rParameters.Has("sub_model_part_name") ? rParameters["sub_model_part_name"].GetString() : ""; if (sub_model_part_name == "") { return rModel.GetModelPart(r_model_part_name); } else { return rModel.GetModelPart(r_model_part_name).GetSubModelPart(sub_model_part_name);}}(rModel, ThisParameters)), mThisParameters(ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)