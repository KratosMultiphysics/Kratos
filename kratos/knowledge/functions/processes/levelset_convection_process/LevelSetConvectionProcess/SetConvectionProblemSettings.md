---
type: C++ Method
title: SetConvectionProblemSettings
resource: processes/levelset_convection_process.h#L507-L529
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/convection_diffusion_settings/SetGradientVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetFillProcessInfoFormulationDataFunction
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetConvectionProblemSettings()`

# Calls

- [SetGradientVariable](../../../../functions/includes/convection_diffusion_settings/SetGradientVariable.md)
- [GetFillProcessInfoFormulationDataFunction](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetFillProcessInfoFormulationDataFunction.md)

# Called by

- [LevelSetConvectionProcess](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess.md)