---
type: C++ Method
title: WriteEmptyElementTag
resource: utilities/xml_utilities/xml_element.cpp#L71-L84
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Write
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_ascii_nd_data_element/XmlAsciiNDDataElement<TDataType>/Write
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_base64_binary_nd_data_element/XmlBase64BinaryNDDataElement<TDataType>/Write
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/Write
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void XmlElement::WriteEmptyElementTag( std::ostream& rOStream, const IndexType Level) const`

# Called by

- [Write](../../../../../functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Write.md)
- [Write](../../../../../functions/utilities/xml_utilities/xml_ascii_nd_data_element/XmlAsciiNDDataElement<TDataType>/Write.md)
- [Write](../../../../../functions/utilities/xml_utilities/xml_base64_binary_nd_data_element/XmlBase64BinaryNDDataElement<TDataType>/Write.md)
- [Write](../../../../../functions/utilities/xml_utilities/xml_elements_array/XmlElementsArray/Write.md)