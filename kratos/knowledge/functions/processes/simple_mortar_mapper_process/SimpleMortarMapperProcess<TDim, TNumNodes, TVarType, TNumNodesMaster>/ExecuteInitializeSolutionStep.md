---
type: C++ Method
title: ExecuteInitializeSolutionStep
resource: processes/simple_mortar_mapper_process.cpp#L164-L196
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/UpdateInterface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>:: ExecuteInitializeSolutionStep()`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [UpdateInterface](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/UpdateInterface.md)
- [ExecuteExplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping.md)
- [ExecuteImplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping.md)