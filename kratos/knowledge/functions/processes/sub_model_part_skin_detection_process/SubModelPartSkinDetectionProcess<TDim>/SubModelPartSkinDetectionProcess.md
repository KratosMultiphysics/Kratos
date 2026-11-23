---
type: C++ Method
title: SubModelPartSkinDetectionProcess
resource: processes/sub_model_part_skin_detection_process.cpp#L57-L86
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SubModelPartSkinDetectionProcess<TDim>::SubModelPartSkinDetectionProcess( ModelPart& rModelPart, Parameters Settings) : SkinDetectionProcess<TDim>(rModelPart, Settings, this->GetDefaultParameters())`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)