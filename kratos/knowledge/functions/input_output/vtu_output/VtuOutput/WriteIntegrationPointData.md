---
type: C++ Method
title: WriteIntegrationPointData
resource: input_output/vtu_output.cpp#L1192-L1365
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
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/WriteData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<std::string, std::string> VtuOutput::WriteIntegrationPointData( TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor, TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor, UnstructuredGridData& rUnstructuredGridData, const std::string& rOutputFileNamePrefix, const IndexType Step) const`

# Calls

- [AddAttribute](../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)
- [GetEndianness](../../../../functions/input_output/vtu_output/GetEndianness.md)
- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [WriteData](../../../../functions/input_output/vtu_output/VtuOutput/WriteData.md)