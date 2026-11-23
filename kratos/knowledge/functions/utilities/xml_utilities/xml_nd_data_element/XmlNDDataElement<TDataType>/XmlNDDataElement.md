---
type: C++ Method
title: XmlNDDataElement
resource: utilities/xml_utilities/xml_nd_data_element.cpp#L29-L61
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`XmlNDDataElement<TDataType>::XmlNDDataElement( const std::string& rDataArrayName, typename NDData<TDataType>::Pointer pNDData) : BaseType("DataArray"), mpNDData(pNDData)`

# Calls

- [AddAttribute](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)