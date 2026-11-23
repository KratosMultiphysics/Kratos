---
type: C++ Method
title: ComputeResidualMatrix
resource: processes/simple_mortar_mapper_process.cpp#L632-L660
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/PerformMortarOperations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ComputeResidualMatrix( Matrix& rResidualMatrix, const GeometryType& rSlaveGeometry, const GeometryType& rMasterGeometry, const MortarOperatorType& rThisMortarOperators )`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [PerformMortarOperations](../../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)