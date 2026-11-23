---
type: C++ Function
title: LengthCheck
resource: utilities/mortar_utilities.cpp#L30-L40
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool LengthCheck( const GeometryPointType& rGeometryLine, const double Tolerance )`

# Called by

- [AssemblyMortarOperators](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [CalculateAe](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/CalculateAe.md)