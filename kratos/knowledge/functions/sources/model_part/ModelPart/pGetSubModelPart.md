---
type: C++ Method
title: pGetSubModelPart
resource: sources/model_part.cpp#L2065-L2080
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart* ModelPart::pGetSubModelPart(std::string const& SubModelPartName)`

# Calls

- [ErrorNonExistingSubModelPart](../../../../functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart.md)

# Called by

- [SensitivityBuilder](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)