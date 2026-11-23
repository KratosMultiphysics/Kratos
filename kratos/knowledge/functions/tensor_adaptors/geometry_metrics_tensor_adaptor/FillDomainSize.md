---
type: C++ Function
title: FillDomainSize
resource: tensor_adaptors/geometry_metrics_tensor_adaptor.cpp#L52-L59
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tensor_adaptors/geometry_metrics_tensor_adaptor/FillData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FillDomainSize( Kratos::span<double> DataSpan, const TContainerType& rContainer)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [FillData](../../../functions/tensor_adaptors/geometry_metrics_tensor_adaptor/FillData.md)