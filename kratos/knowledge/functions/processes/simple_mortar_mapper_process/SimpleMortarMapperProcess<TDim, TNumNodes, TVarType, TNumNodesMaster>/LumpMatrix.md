---
type: C++ Method
title: LumpMatrix
resource: processes/simple_mortar_mapper_process.cpp#L550-L560
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/PerformMortarOperations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::LumpMatrix(BoundedMatrixType& rInputMatrix)`

# Called by

- [PerformMortarOperations](../../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)