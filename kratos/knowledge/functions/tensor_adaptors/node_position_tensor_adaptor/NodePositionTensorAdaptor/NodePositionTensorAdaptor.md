---
type: C++ Method
title: NodePositionTensorAdaptor
resource: tensor_adaptors/node_position_tensor_adaptor.cpp#L26-L49
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`NodePositionTensorAdaptor::NodePositionTensorAdaptor( ModelPart::NodesContainerType::Pointer pContainer, Globals::Configuration Configuration, const std::vector<unsigned int>& rDataShape) : mConfiguration(Configuration)`

# Calls

- [DataShape](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/DataShape.md)