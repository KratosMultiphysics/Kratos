---
type: C++ Method
title: CollectShapeFunctions
resource: tensor_adaptors/geometries_tensor_adaptor.cpp#L141-L166
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/FillData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometriesTensorAdaptor::CollectShapeFunctions( double* pData, const TContainerType& rContainer, GeometryData::IntegrationMethod Method, const DenseVector<unsigned int>& Shape)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [FillData](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/FillData.md)