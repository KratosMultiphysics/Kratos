---
type: C++ Function
title: SetGradientVariable
resource: includes/convection_diffusion_settings.h#L207-L211
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/levelset_convection_process/ProcessInfoDataContainer/RestoreProcessInfoData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetGradientVariable(const Variable<array_1d<double,3> >& rvar)`

# Called by

- [RestoreProcessInfoData](../../../functions/processes/levelset_convection_process/ProcessInfoDataContainer/RestoreProcessInfoData.md)
- [Execute](../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)
- [SetConvectionProblemSettings](../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings.md)