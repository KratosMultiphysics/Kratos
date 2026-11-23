---
type: C++ Method
title: GetIntegrationMethod
resource: tensor_adaptors/geometries_tensor_adaptor.cpp#L30-L45
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometryData::IntegrationMethod GeometriesTensorAdaptor::GetIntegrationMethod( const TContainerType& rContainer, GeometryData::IntegrationMethod UserMethod)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)