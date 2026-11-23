---
type: C++ Method
title: WriteUnstructuredGridData
resource: input_output/vtu_output.cpp#L1087-L1189
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetEndianness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddConnectivityData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddFields
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<std::string, std::string> VtuOutput::WriteUnstructuredGridData( TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor, TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor, UnstructuredGridData& rUnstructuredGridData, const std::string& rOutputFileNamePrefix, const IndexType Step) const`

# Calls

- [AddAttribute](../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)
- [GetEndianness](../../../../functions/input_output/vtu_output/GetEndianness.md)
- [pGetStorage](../../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/pGetStorage.md)
- [AddConnectivityData](../../../../functions/input_output/vtu_output/AddConnectivityData.md)
- [AddFields](../../../../functions/input_output/vtu_output/AddFields.md)
- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)