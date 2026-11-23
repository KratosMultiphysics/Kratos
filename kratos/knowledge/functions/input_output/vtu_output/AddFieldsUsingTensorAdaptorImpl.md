---
type: C++ Function
title: AddFieldsUsingTensorAdaptorImpl
resource: input_output/vtu_output.cpp#L340-L433
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/GetWritingNDData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddFieldsUsingTensorAdaptorImpl( XmlElementsArray& rXmlElement, TContainerPointerType pContainer, const Variable<TDataType>& rVariable, TXmlDataElementWrapper& rXmlDataElementWrapper, const DataCommunicator& rDataCommunicator, const std::vector<IndexType>& rWritingIndices, const IndexType EchoLevel, TArgs&&... rArgs)`

# Calls

- [GetWritingNDData](../../../functions/input_output/vtu_output/GetWritingNDData.md)
- [pGetStorage](../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage.md)