---
type: C++ Method
title: CalculateEnrichedShapeFuncionsExtended
resource: utilities/enrichment_utilities.h#L1230-L1634
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/active_node_in_replacement_shape_function
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int CalculateEnrichedShapeFuncionsExtended(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (4)>& NEnriched)`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [active_node_in_replacement_shape_function](../../../../functions/utilities/enrichment_utilities/active_node_in_replacement_shape_function.md)