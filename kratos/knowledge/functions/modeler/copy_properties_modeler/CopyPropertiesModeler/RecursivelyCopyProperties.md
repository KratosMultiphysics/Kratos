---
type: C++ Method
title: RecursivelyCopyProperties
resource: modeler/copy_properties_modeler.cpp#L79-L92
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CopyPropertiesModeler::RecursivelyCopyProperties( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)`

# Calls

- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/SetupModelPart.md)