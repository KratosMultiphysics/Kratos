---
type: C++ Method
title: GetShape
resource: tensor_adaptors/geometries_tensor_adaptor.cpp#L48-L119
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseVector<unsigned int> GeometriesTensorAdaptor::GetShape( const TContainerType& rContainer, GeometriesTensorAdaptor::DatumType Datum, GeometryData::IntegrationMethod Method)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)