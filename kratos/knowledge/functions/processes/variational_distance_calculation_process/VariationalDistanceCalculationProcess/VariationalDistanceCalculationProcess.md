---
type: C++ Method
title: VariationalDistanceCalculationProcess
resource: processes/variational_distance_calculation_process.h#L132-L142
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/AsFalse
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/InitializeSolutionStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`VariationalDistanceCalculationProcess( Model& rModel, typename TLinearSolver::Pointer pLinearSolver, Parameters ThisParameters) : VariationalDistanceCalculationProcess( rModel, pLinearSolver, Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver), ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [AsFalse](../../../../functions/containers/flags/AsFalse.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [ReGenerateDistanceModelPart](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart.md)
- [InitializeSolutionStrategy](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/InitializeSolutionStrategy.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)