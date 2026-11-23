---
type: C++ Method
title: CalculateDiscontinuousShapeFunctions
resource: utilities/discont_utils.h#L435-L672
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateDiscontinuousShapeFunctions(const BoundedMatrix<double,(2+1), 2 >& rPoints, BoundedMatrix<double, (2+1), 2 >& rDN_DX, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2+1)>& rNenriched, array_1d<double,3>& rEdgeAreas)`