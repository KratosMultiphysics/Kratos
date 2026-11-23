---
type: C++ Method
title: XmlAsciiNDDataElement
resource: utilities/xml_utilities/xml_ascii_nd_data_element.cpp#L29-L38
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

`XmlAsciiNDDataElement<TDataType>::XmlAsciiNDDataElement( const std::string& rDataArrayName, typename NDData<TDataType>::Pointer pNDData, const IndexType Precision) : BaseType(rDataArrayName, pNDData), mPrecision(Precision)`

# Calls

- [AddAttribute](../../../../../functions/utilities/xml_utilities/xml_element/XmlElement/AddAttribute.md)