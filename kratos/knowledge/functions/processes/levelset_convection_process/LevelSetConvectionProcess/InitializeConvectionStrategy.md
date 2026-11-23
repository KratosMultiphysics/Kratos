---
type: C++ Method
title: InitializeConvectionStrategy
resource: processes/levelset_convection_process.h#L1037-L1076
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/ReGenerateConvectionModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeConvectionStrategy(BuilderAndSolverPointerType pBuilderAndSolver)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ReGenerateConvectionModelPart](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/ReGenerateConvectionModelPart.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [LevelSetConvectionProcess](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess.md)