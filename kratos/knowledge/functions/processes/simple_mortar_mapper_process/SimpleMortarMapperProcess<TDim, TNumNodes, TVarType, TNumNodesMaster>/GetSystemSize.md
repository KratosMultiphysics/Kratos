---
type: C++ Method
title: GetSystemSize
resource: processes/simple_mortar_mapper_process.cpp#L566-L569
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetSystemSize(std::size_t& rSizeSystem)`

# Called by

- [ExecuteExplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping.md)