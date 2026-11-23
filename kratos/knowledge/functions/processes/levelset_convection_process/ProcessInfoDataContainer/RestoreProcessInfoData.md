---
type: C++ Method
title: RestoreProcessInfoData
resource: processes/levelset_convection_process.h#L85-L91
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/convection_diffusion_settings/SetGradientVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RestoreProcessInfoData(ProcessInfo& rOutputProcessInfo) const`

# Calls

- [SetGradientVariable](../../../../functions/includes/convection_diffusion_settings/SetGradientVariable.md)

# Called by

- [Execute](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)