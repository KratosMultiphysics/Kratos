---
type: C++ Method
title: Write
resource: utilities/xml_utilities/xml_appended_data_element_wrapper.cpp#L114-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/WriteEmptyElementTag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/WriteElementTagStart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/WriteElementTagEnd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void XmlAppendedDataElementWrapper::Write( std::ostream& rOStream, const IndexType Level) const`

# Calls

- [WriteEmptyElementTag](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/WriteEmptyElementTag.md)
- [WriteElementTagStart](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/WriteElementTagStart.md)
- [WriteElementTagEnd](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/WriteElementTagEnd.md)