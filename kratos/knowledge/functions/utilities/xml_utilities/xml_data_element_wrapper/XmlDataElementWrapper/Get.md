---
type: C++ Method
title: Get
resource: utilities/xml_utilities/xml_data_element_wrapper.cpp#L30-L55
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_utils/XmlUtilities/AddDataArrayAttributes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`XmlElement::Pointer XmlDataElementWrapper::Get( const std::string& rDataArrayName, NDDataPointerType pNDData)`

# Calls

- [AddAttribute](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)
- [visit](../../../../../functions/python/variable_indexing_python/visit.md)
- [AddDataArrayAttributes](../../../../../functions/utilities/xml_utilities/xml_utils/XmlUtilities/AddDataArrayAttributes.md)