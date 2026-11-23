---
type: C++ Function
title: CalculateTriangleDiscontinuousShapeFunctions_ZeroInBoundary
resource: utilities/discont_2d.h#L312-L657
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateTriangleDiscontinuousShapeFunctions_ZeroInBoundary(TMatrixType const& rPoints, TGradientType const& DN_DX, TVectorType rDistances, TVectorType& rVolumes, TMatrixType& rGPShapeFunctionValues, TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched, //and information about the interfase: TVectorType& face_gauss_N, TVectorType& face_gauss_Nenriched, double& face_Area, TVectorType& face_n ,unsigned int& type_of_cut)`