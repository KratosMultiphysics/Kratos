---
type: C++ Function
title: GaussPointUnitNormal
resource: utilities/mortar_utilities.cpp#L113-L128
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double,3> GaussPointUnitNormal( const Vector& rN, const GeometryType& rGeometry )`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [AssemblyMortarOperators](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)