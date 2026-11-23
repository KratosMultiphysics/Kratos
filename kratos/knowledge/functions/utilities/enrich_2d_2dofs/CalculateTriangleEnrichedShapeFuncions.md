---
type: C++ Function
title: CalculateTriangleEnrichedShapeFuncions
resource: utilities/enrich_2d_2dofs.h#L86-L407
generated:
  by: okf-rs/0.7.0
---

# Signature

`static int CalculateTriangleEnrichedShapeFuncions(TMatrixType const& rPoints, TGradientType const& DN_DX, TVectorType rDistances, TVectorType& rVolumes, TMatrixType& rGPShapeFunctionValues, TVectorType& rPartitionsSign, std::vector<TMatrixType>& rGradientsValue, TMatrixType& NEnriched, //and information about the interfase: TVectorType& face_gauss_N, TVectorType& face_gauss_Nenriched, double& face_Area, TVectorType& face_n ,unsigned int& type_of_cut)`