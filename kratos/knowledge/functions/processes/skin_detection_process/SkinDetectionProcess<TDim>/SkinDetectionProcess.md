---
type: C++ Method
title: SkinDetectionProcess
resource: processes/skin_detection_process.cpp#L27-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SkinDetectionProcess<TDim>::SkinDetectionProcess( ModelPart& rModelPart, Parameters ThisParameters ) : mrModelPart(rModelPart), mThisParameters(ThisParameters)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)