---
type: C++ Method
title: CollectFixityData
resource: tensor_adaptors/fixity_tensor_adaptor.cpp#L112-L124
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FixityTensorAdaptor::CollectFixityData( Kratos::span<bool> Span, const std::vector<const Variable<double>*>& rDofsVarPointerList, const ModelPart::NodesContainerType& rNodes)`

# Called by

- [CollectData](../../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CollectData.md)