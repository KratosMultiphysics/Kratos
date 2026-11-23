---
type: C++ Method
title: CreateSlaveConectivityDatabase
resource: processes/simple_mortar_mapper_process.cpp#L575-L592
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CreateSlaveConectivityDatabase( std::size_t& rSizeSystem, IntMap& rConectivityDatabase, IntMap& rInverseConectivityDatabase )`

# Called by

- [ExecuteImplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping.md)