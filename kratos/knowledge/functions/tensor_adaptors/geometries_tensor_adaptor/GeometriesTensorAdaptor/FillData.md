---
type: C++ Method
title: FillData
resource: tensor_adaptors/geometries_tensor_adaptor.cpp#L249-L270
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctionsDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectJacobians
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometriesTensorAdaptor::FillData(double* pData, const TContainerType& rContainer, GeometriesTensorAdaptor::DatumType Datum, GeometryData::IntegrationMethod Method, const DenseVector<unsigned int>& Shape)`

# Calls

- [CollectShapeFunctions](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctions.md)
- [CollectShapeFunctionsDerivatives](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctionsDerivatives.md)
- [CollectNormals](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectNormals.md)
- [CollectJacobians](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectJacobians.md)
- [CollectIntegrationWeights](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights.md)