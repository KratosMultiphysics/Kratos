---
type: Python Method
title: execute_test_model_part_io_read_model_part
resource: tests/test_model_part_io.py#L60-L220
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfTables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def execute_test_model_part_io_read_model_part(self, input_mdpa):`

# Calls

- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [NumberOfTables](../../../../functions/includes/model_part/NumberOfTables.md)
- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [test_model_part_io_read_model_part](../../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part.md)