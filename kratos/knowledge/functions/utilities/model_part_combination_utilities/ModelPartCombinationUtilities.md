---
type: C++ Function
title: ModelPartCombinationUtilities
resource: utilities/model_part_combination_utilities.h#L63-L66
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts_fail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartCombinationUtilities(Model& rModel): mrModel(rModel)`

# Called by

- [SetupModelPart](../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [test_combine_model_parts](../../../functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts.md)
- [test_combine_model_parts_fail](../../../functions/tests/test_model_part_combination_utilities/TestModelPartCombinationUtilities/test_combine_model_parts_fail.md)