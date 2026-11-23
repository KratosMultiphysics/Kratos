---
type: C++ Method
title: CollectIds
resource: tensor_adaptors/connectivity_ids_tensor_adaptor.cpp#L64-L81
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConnectivityIdsTensorAdaptor::CollectIds(Kratos::span<int> Span, const TContainerType& rContainer, const std::size_t NumNodes)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CollectData](../../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectData.md)