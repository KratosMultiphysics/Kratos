---
type: C++ Method
title: LevelSetConvectionProcess
resource: processes/levelset_convection_process.h#L130-L139
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`LevelSetConvectionProcess( Model& rModel, typename TLinearSolver::Pointer pLinearSolver, Parameters ThisParameters) : LevelSetConvectionProcess( rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), pLinearSolver, ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [InitializeConvectionStrategy](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetConvectionElementDefaultParameters](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementDefaultParameters.md)
- [SetConvectionProblemSettings](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/SetConvectionProblemSettings.md)