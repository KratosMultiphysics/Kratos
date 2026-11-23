---
type: Python Function
title: ReadModelPart
resource: python_scripts/testing/utilities.py#L22-L39
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/utilities/ReadDistributedModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/ReadSerialModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ReadModelPart(mdpa_file_name, model_part, settings=None):`

# Calls

- [ReadDistributedModelPart](../../../../functions/python_scripts/testing/utilities/ReadDistributedModelPart.md)
- [ReadSerialModelPart](../../../../functions/python_scripts/testing/utilities/ReadSerialModelPart.md)