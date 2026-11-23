---
type: C++ Method
title: CalculateUnitNormal
resource: containers/edge_based_data_structure.h#L588-L612
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateUnitNormal( const Geometry<Node>& rGeometry, array_1d<double, TDim>& rUnitNormal) const`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)