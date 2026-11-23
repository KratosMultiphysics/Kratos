---
type: C++ Method
title: ExecuteExplicitMapping
resource: processes/simple_mortar_mapper_process.cpp#L726-L892
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetSystemSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ResetNodalArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetReferenceArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckAndPerformSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CreateInverseDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckWholeVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ExecuteExplicitMapping()`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetSystemSize](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetSystemSize.md)
- [ResetNodalArea](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ResetNodalArea.md)
- [GetReferenceArea](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetReferenceArea.md)
- [CheckAndPerformSearch](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckAndPerformSearch.md)
- [CreateInverseDatabase](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CreateInverseDatabase.md)
- [CheckWholeVector](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CheckWholeVector.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteInitializeSolutionStep.md)