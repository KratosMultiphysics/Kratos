---
type: C++ Function
title: CopyPropertiesValues
resource: utilities/properties_utilities.cpp#L24-L39
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_properties/TestProperties/test_properties_utilities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CopyPropertiesValues( const Properties& rOriginProperties, Properties& rDestinatiionProperties )`

# Calls

- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [test_properties_utilities](../../../functions/tests/test_properties/TestProperties/test_properties_utilities.md)