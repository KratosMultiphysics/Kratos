---
type: C++ Function
title: AddFields
resource: input_output/vtu_output.cpp#L482-L561
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetWritingNDData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddFields( XmlElementsArray& rXmlElement, const std::map<std::string, VtuOutput::SupportedTensorAdaptorPointerType>& rMap, TXmlDataElementWrapper& rXmlDataElementWrapper, const std::vector<IndexType>& rWritingIndices, const DataCommunicator& rDataCommunicator, const IndexType EchoLevel)`

# Calls

- [visit](../../../functions/python/variable_indexing_python/visit.md)
- [GetWritingNDData](../../../functions/input_output/vtu_output/GetWritingNDData.md)
- [pGetStorage](../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage.md)

# Called by

- [WriteUnstructuredGridData](../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)