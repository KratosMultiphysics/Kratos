---
type: C++ Method
title: TimeAveragingProcess
resource: processes/time_averaging_process.cpp#L28-L66
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/process_factory/Factory
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TimeAveragingProcess::TimeAveragingProcess(Model& rModel, Parameters rParameters) : mrModel(rModel), mrParameters(rParameters)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [Factory](../../../../functions/python_scripts/process_factory/Factory.md)