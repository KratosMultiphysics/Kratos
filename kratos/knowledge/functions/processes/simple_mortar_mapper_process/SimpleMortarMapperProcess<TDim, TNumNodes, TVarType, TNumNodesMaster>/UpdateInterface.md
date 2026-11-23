---
type: C++ Method
title: UpdateInterface
resource: processes/simple_mortar_mapper_process.cpp#L1145-L1170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/UpdateInterfaceWithoutArgs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/UpdateInterfaceWithOptions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mapper_to_python/UpdateInterfaceWithSearchRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::UpdateInterface()`

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteInitializeSolutionStep.md)
- [UpdateInterfaceWithoutArgs](../../../../functions/python/add_mapper_to_python/UpdateInterfaceWithoutArgs.md)
- [UpdateInterfaceWithOptions](../../../../functions/python/add_mapper_to_python/UpdateInterfaceWithOptions.md)
- [UpdateInterfaceWithSearchRadius](../../../../functions/python/add_mapper_to_python/UpdateInterfaceWithSearchRadius.md)