---
type: C++ Method
title: AddAttribute
resource: utilities/xml_utilities/xml_element.cpp#L32-L37
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_xml_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/XmlAppendedDataElementWrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Get
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_ascii_nd_data_element/XmlAsciiNDDataElement<TDataType>/XmlAsciiNDDataElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_base64_binary_nd_data_element/XmlBase64BinaryNDDataElement<TDataType>/XmlBase64BinaryNDDataElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_data_element_wrapper/XmlDataElementWrapper/Get
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_nd_data_element/XmlNDDataElement<TDataType>/XmlNDDataElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_utils/constexpr
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void XmlElement::AddAttribute( const std::string& rName, const std::string& rValue)`

# Called by

- [WriteUnstructuredGridData](../../../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [WriteIntegrationPointData](../../../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [PrintOutput](../../../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_xml_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [XmlAppendedDataElementWrapper](../../../../../functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/XmlAppendedDataElementWrapper.md)
- [Get](../../../../../functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Get.md)
- [XmlAsciiNDDataElement](../../../../../functions/utilities/xml_utilities/xml_ascii_nd_data_element/XmlAsciiNDDataElement<TDataType>/XmlAsciiNDDataElement.md)
- [XmlBase64BinaryNDDataElement](../../../../../functions/utilities/xml_utilities/xml_base64_binary_nd_data_element/XmlBase64BinaryNDDataElement<TDataType>/XmlBase64BinaryNDDataElement.md)
- [Get](../../../../../functions/utilities/xml_utilities/xml_data_element_wrapper/XmlDataElementWrapper/Get.md)
- [XmlNDDataElement](../../../../../functions/utilities/xml_utilities/xml_nd_data_element/XmlNDDataElement<TDataType>/XmlNDDataElement.md)
- [constexpr](../../../../../functions/utilities/xml_utilities/xml_utils/constexpr.md)