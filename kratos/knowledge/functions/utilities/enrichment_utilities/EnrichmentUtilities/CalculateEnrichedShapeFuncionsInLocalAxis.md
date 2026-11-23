---
type: C++ Method
title: CalculateEnrichedShapeFuncionsInLocalAxis
resource: utilities/enrichment_utilities.h#L1996-L2357
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int CalculateEnrichedShapeFuncionsInLocalAxis(BoundedMatrix<double,(2+1), 2 >& rOriginalPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX_original, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2)>& NEnriched, BoundedMatrix<double,(2), 2 >& rRotationMatrix, BoundedMatrix<double, (2+1), 2 > & DN_DX_in_local_axis)`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)