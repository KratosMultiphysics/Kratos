---
type: C++ Method
title: CalculateEnrichedShapeFuncions
resource: utilities/enrichment_utilities.h#L1638-L1991
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateEnrichedShapeFuncions(BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2)>& NEnriched, //information about the interfase: Standard shape functions values, Enrichment Functions (first(gradient), second (positive_side,negative side) ) and Area array_1d<double,(3)>& rGPShapeFunctionValues_in_interfase, array_1d<double,(3)>& NEnriched_in_interfase, double & InterfaseArea)`