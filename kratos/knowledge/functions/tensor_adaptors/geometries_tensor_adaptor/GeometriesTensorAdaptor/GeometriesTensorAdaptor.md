---
type: C++ Method
title: GeometriesTensorAdaptor
resource: tensor_adaptors/geometries_tensor_adaptor.cpp#L272-L297
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometriesTensorAdaptor::GeometriesTensorAdaptor( ContainerPointerType pContainer, DatumType Datum, GeometryData::IntegrationMethod IntegrationMethod) : TensorAdaptor<double>(), mDatum(Datum), mIntegrationMethod(IntegrationMethod)`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)