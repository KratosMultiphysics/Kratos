---
type: C++ Method
title: ComputeNodalGradientProcess
resource: processes/compute_nodal_gradient_process.cpp#L29-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CheckVariableExists
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ComputeNodalGradientProcess<TOutputHistorical>::ComputeNodalGradientProcess( Model& rModel, Parameters ThisParameters ) : mrModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CheckVariableExists](../../../../functions/utilities/variable_utils/CheckVariableExists.md)