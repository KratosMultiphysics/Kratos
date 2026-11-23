---
type: C++ Method
title: CalculateEnrichedShapeFuncions
resource: utilities/enrichment_utilities.h#L911-L1227
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateEnrichedShapeFuncions( BoundedMatrix<double,3, 2 >& rPoints, BoundedMatrix<double, 3, 2 >& DN_DX, array_1d<double,3>& rDistances, array_1d<double,3>& rVolumes, BoundedMatrix<double, 3, 3 >& rGPShapeFunctionValues, array_1d<double,3>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3, 2>& NEnriched)`