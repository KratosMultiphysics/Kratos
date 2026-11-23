---
type: C++ Method
title: StoreFixityData
resource: tensor_adaptors/fixity_tensor_adaptor.cpp#L126-L142
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FixityTensorAdaptor::StoreFixityData( Kratos::span<bool> Span, const std::vector<const Variable<double>*>& rDofsVarPointerList, const ModelPart::NodesContainerType& rNodes)`

# Calls

- [Fix](../../../../functions/includes/node/Fix.md)
- [Free](../../../../functions/includes/node/Free.md)

# Called by

- [StoreData](../../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreData.md)