---
type: C++ Method
title: AddConvectiveBoundaryValue
resource: containers/edge_based_data_structure.h#L143-L156
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddConvectiveBoundaryValue( const double Weight, const double Ni, const double Nj, const array_1d<double,TDim>& rUnitNormal)`

# Called by

- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)