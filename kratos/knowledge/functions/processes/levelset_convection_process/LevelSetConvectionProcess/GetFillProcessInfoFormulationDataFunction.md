---
type: C++ Method
title: GetFillProcessInfoFormulationDataFunction
resource: processes/levelset_convection_process.h#L1014-L1035
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const virtual std::function<void(ModelPart&)> GetFillProcessInfoFormulationDataFunction()`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)

# Called by

- [Execute](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)
- [SetConvectionProblemSettings](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings.md)