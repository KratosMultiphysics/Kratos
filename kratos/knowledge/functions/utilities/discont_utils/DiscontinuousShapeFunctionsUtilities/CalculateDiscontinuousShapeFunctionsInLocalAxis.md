---
type: C++ Method
title: CalculateDiscontinuousShapeFunctionsInLocalAxis
resource: utilities/discont_utils.h#L675-L927
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateDiscontinuousShapeFunctionsInLocalAxis( BoundedMatrix<double,(2+1), 2 >& rOriginalPoints, BoundedMatrix<double,(2+1), 2 >& rRotatedPoints, BoundedMatrix<double, (2+1), 2 >& DN_DX_original, array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, BoundedMatrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, BoundedMatrix<double,3*(2-1), (2+1)>& Nenriched, BoundedMatrix<double,(2), 2 >& rRotationMatrix, BoundedMatrix<double, (2+1), 2 >& DN_DX_in_local_axis )`