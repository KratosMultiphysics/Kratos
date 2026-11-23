---
type: C++ Method
title: ReGenerateConvectionModelPart
resource: processes/levelset_convection_process.h#L531-L600
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeDistanceModelPartDatabases
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart)`

# Calls

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [InitializeDistanceModelPartDatabases](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeDistanceModelPartDatabases.md)

# Called by

- [Execute](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)
- [InitializeConvectionStrategy](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy.md)