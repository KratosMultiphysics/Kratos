---
type: C++ Method
title: CopyCommonData
resource: modeler/connectivity_preserve_modeler.cpp#L224-L257
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/pGetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConnectivityPreserveModeler::CopyCommonData( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [pGetNodalSolutionStepVariablesList](../../../../functions/includes/model_part/pGetNodalSolutionStepVariablesList.md)