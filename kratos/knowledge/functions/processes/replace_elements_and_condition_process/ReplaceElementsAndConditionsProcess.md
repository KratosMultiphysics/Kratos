---
type: C++ Function
title: ReplaceElementsAndConditionsProcess
resource: processes/replace_elements_and_condition_process.h#L60-L73
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/WriteJsonString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ReplaceElementsAndConditionsProcess( Model& rModel, Parameters Settings ) : Process(Flags()) , mrModelPart(rModel.GetModelPart(Settings["model_part_name"].GetString())), mSettings( Settings.WriteJsonString())`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [WriteJsonString](../../../functions/sources/kratos_parameters/Parameters/WriteJsonString.md)