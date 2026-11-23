---
type: C++ Method
title: SimpleMortarMapperProcess
resource: processes/simple_mortar_mapper_process.cpp#L46-L74
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters ThisParameters, LinearSolverType::Pointer pThisLinearSolver ): mOriginModelPart(rOriginModelPart), mDestinationModelPart(rDestinationModelPart), mpOriginVariable(&(KratosComponents<TVarType>::Get(ThisParameters["origin_variable"].GetString()))), mpDestinationVariable((ThisParameters["destination_variable"].GetString() == "") ? mpOriginVariable : &(KratosComponents<TVarType>::Get(ThisParameters["destination_variable"].GetString()))), mThisParameters(ThisParameters), mpThisLinearSolver(pThisLinearSolver)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)