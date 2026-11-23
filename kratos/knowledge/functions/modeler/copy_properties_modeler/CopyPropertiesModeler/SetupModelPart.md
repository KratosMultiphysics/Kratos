---
type: C++ Method
title: SetupModelPart
resource: modeler/copy_properties_modeler.cpp#L61-L77
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/RecursivelyCopyProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/ReplaceProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CopyPropertiesModeler::SetupModelPart()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyCopyProperties](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/RecursivelyCopyProperties.md)
- [ReplaceProperties](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/ReplaceProperties.md)