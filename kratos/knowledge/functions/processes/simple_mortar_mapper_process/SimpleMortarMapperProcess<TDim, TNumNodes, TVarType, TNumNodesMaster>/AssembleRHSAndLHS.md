---
type: C++ Method
title: AssembleRHSAndLHS
resource: processes/simple_mortar_mapper_process.cpp#L666-L697
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

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AssembleRHSAndLHS( MatrixType& rA, std::vector<VectorType>& rb, const SizeType VariableSize, const Matrix& rResidualMatrix, const GeometryType& rSlaveGeometry, IntMap& rInverseConectivityDatabase, const MortarOperatorType& rThisMortarOperators )`

# Called by

- [PerformMortarOperations](../../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)