---
type: C++ Method
title: ProcessInfoDataContainer
resource: processes/levelset_convection_process.h#L78-L83
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ProcessInfoDataContainer(const ProcessInfo& rInputProcessInfo) : DeltaTime(rInputProcessInfo.GetValue(DELTA_TIME)) , pUnknownVariable(&((rInputProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS))->GetUnknownVariable())) , pGradientVariable(&((rInputProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS))->GetGradientVariable())) , pConvectionVariable(&((rInputProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS))->GetConvectionVariable()))`

# Called by

- [Execute](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)